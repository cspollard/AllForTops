{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}

module Main where

import           Atlas
import           Atlas.CrossSections
import           Atlas.Histogramming
import           Control.Applicative           (ZipList (..))
import           Control.Comonad               (duplicate, extract)
import qualified Control.Foldl                 as F
import           Control.Lens                  (over, view)
import           Control.Monad                 (forM_)
import           Control.Monad.Reader
import           Control.Monad.Trans.Maybe
import           Control.Monad.Writer
import           Control.Monad.Writer.Class
import           Data.Functor.Compose
import           Data.Histogram.Bin.Transform
import qualified Data.Histogram.Generic        as H
import qualified Data.IntMap                   as IM
import           Data.Maybe                    (fromMaybe)
import           Data.Monoid                   ((<>))
import qualified Data.Text                     as T
import           Data.TFile
import           Data.TTree
import qualified Data.Vector                   as V
import           Debug.Trace
import           GHC.Float
import           Options.Generic
import           Pipes
import           Pipes.Lift
import qualified Pipes.Prelude                 as P
import           System.IO
import           System.Random.MWC.Probability

data Args =
  Args
  { outfolder :: String
  , infiles   :: String
  , xsecfile  :: String
  , lumi      :: Double
  } deriving (Show, Generic)

instance ParseRecord Args where


data LargeJet =
  LargeJet
  { ljFourMom :: PtEtaPhiE
  , ljMass    :: Double
  } deriving Show



readFourMoms
  :: (MonadIO m, MonadThrow m)
  => String -> TreeRead m (ZipList PtEtaPhiE)
readFourMoms prefix = do
  pts <- fmap ((/1e3) . float2Double) <$> readBranch (prefix ++ "pt")
  etas <- fmap float2Double <$> readBranch (prefix ++ "eta")
  phis <- fmap float2Double <$> readBranch (prefix ++ "phi")
  es <- fmap ((/1e3) . float2Double) <$> readBranch (prefix ++ "e")

  return $ PtEtaPhiE <$> pts <*> etas <*> phis <*> es


readLargeJets :: (MonadIO m, MonadThrow m) => TreeRead m [LargeJet]
readLargeJets = do
  ms <- fmap ((/1e3) . float2Double) <$> readBranch "ljet_m"
  ljFourMoms <- readFourMoms "ljet_"

  return . getZipList $ LargeJet <$> ljFourMoms <*> ms


newtype RCJet =
  RCJet
  { rcFourMom :: PtEtaPhiE
  } deriving Show


readRCJets :: (MonadIO m, MonadThrow m) => TreeRead m [RCJet]
readRCJets = (fmap RCJet . getZipList) <$> readFourMoms "rcjet_"


data Jet =
  Jet
  { jFourMom  :: PtEtaPhiE
  , jMV2c1077 :: Bool
  } deriving Show


readJets :: (MonadIO m, MonadThrow m) => TreeRead m [Jet]
readJets = do
  fourMoms <- readFourMoms "jet_"
  btags <- fmap cToB <$> readBranch "jet_isbtagged_MV2c10_77"

  return . getZipList $ Jet <$> fourMoms <*> btags

  where
    cToB :: CChar -> Bool
    cToB 0 = False
    cToB _ = True


data Event =
  Event
  { eWeight         :: Double
  , eInclusiveJets  :: [Jet]
  , eAdditionalJets :: [Jet]
  , eTopJets        :: [LargeJet]
  , eLeptons        :: [PtEtaPhiE]
  , eMET            :: PtEtaPhiE
  } deriving Show


readEvent
  :: (MonadIO m, MonadThrow m)
  => Bool -> TreeRead m (Maybe Event)
readEvent isData = do
  pass <-
    or <$> traverse (fmap iToB . readBranch)
      [ "boosted_ejets_2015_MV2c10"
      , "boosted_mujets_2015_MV2c10"
      , "boosted_ejets_2016_MV2c10"
      , "boosted_mujets_2016_MV2c10"
      ]

  if not pass
    then return Nothing
    else fmap Just $ do
      wgt <-
        if isData
          then return 1.0
          else float2Double <$> readBranch "weight_mc"

      jets <- readJets
      topJets <- take 2 . filter topTagged <$> readLargeJets

      let bjets = filter bTagged jets
          tjp4s = ljFourMom <$> topJets
          jets' = removeOverlap tjp4s jets

      return $ Event wgt jets jets' topJets mempty mempty

  where
    removeOverlap ljp4s =
      filter (\j -> all (\lj -> lvDRRap lj (jFourMom j) > 1.0) ljp4s)

    bTagged (Jet _ tagged) = tagged

    topTagged (LargeJet p4 m) = view lvPt p4 > 250 && m > 100

    iToB :: CInt -> Bool
    iToB 0 = False
    iToB _ = True


mJJ :: Hist1DFill LogBinD Double
mJJ = F.premap (,1.0) $ hist1DFill hmJJ
  where
    hmJJ = H.histogramUO (logBinD 600 100 6e3) Nothing (V.replicate 100 mempty)


toHandleF :: MonadIO m => String -> F.FoldM m String ()
toHandleF fn = F.FoldM step start done
  where
    step h s = liftIO (hPutStrLn h s) >> return h
    start =
      liftIO $ do
        h <- openFile fn WriteMode
        hSetBuffering h LineBuffering
        return h

    done h = liftIO $ hClose h


channelF
  :: MonadIO m
  => String
  -> (a -> Maybe b)
  -> F.FoldM m b c
  -> F.FoldM m a (String, c)
channelF s f fol = F.premapM f $ (s,) <$> F.handlesM F.folded fol


sampleEvent :: (MonadIO m, F.PrimMonad m) => Event -> Prob m Bool
sampleEvent Event{..} =
  if eWeight >= 1
    then do
      liftIO $ hPutStrLn stderr "warning: event weight > 1!"
      return True
    else do
      x <- uniform
      return $ eWeight > x


channels :: MonadIO m => String -> F.FoldM m Event [(String, Hist1D LogBinD)]
channels prefix =
  traverse (\(a, b) -> channelF a b (mJJ' a))
  [ (prefix ++ "eq3j_eq2b", eventCut (== 3) (== 2))
  , (prefix ++ "eq3j_eq3b", eventCut (== 3) (== 3))
  , (prefix ++ "eq3j_ge4b", eventCut (== 3) (>= 4))
  , (prefix ++ "ge4j_eq2b", eventCut (>= 4) (== 2))
  , (prefix ++ "ge4j_eq3b", eventCut (>= 4) (== 3))
  , (prefix ++ "ge4j_ge4b", eventCut (>= 4) (>= 4))
  ]

  where
    mJJ' :: MonadIO m => String -> F.FoldM m Double (Hist1D LogBinD)
    mJJ' s =
      const
      <$> F.generalize mJJ
      <*> F.premapM (\m -> "1.0, " ++ show m) (toHandleF s)

    tagged (Jet _ t) = t

    eventCut cutj cutb Event{..} = do
      let nb = foldl (\s j -> if tagged j then s+1 else s) 0 eInclusiveJets
          nj = length eAdditionalJets
          nJ = length eTopJets

      guard $ nJ >= 2
      guard $ cutj nj
      guard $ cutb nb

      let (tj1:tj2:_) = ljFourMom <$> eTopJets

      return . view lvM $ tj1 <> tj2

instance MonadThrow m => MonadThrow (Prob m) where
  throwM = lift . throwM


main :: IO ()
main = do
  args <- getRecord "run-4t" :: IO Args

  hSetBuffering stdout LineBuffering

  (files :: [TFile]) <-
    P.toListM
    $ linesP (infiles args) >-> P.mapM tfileOpen

  (Just (dsid :: CInt)) <- P.head $ for (each files) readDSID

  (hists :: [(String, Hist1D LogBinD)]) <-
    if dsid == 0
      then
        F.impurely P.foldM (channels $ outfolder args ++ "/")
        $ for (each files) (readEvents True)

      else do
        sow <-
          F.impurely P.foldM (F.generalize F.sum)
          $ for (each files) readSumWeights

        putStrLn $ "sum of weights in all files: " ++ show sow

        xsecs <-
          fromMaybe (error "failed to read xsecs") <$> readXSecFile (xsecfile args)

        let (xsec, _) = xsecs IM.! fromEnum dsid
            scale = xsec * lumi args / sow

        hists' <-
          withSystemRandom . asGenIO . sample
          . F.impurely P.foldM (channels $ outfolder args ++ "/")
          $ for (each files) (readEvents False)
            >-> P.map (scaleWgt scale)
            >-> P.filterM sampleEvent

        return
          $ over (traverse.traverse) (scaling scale) hists'

  withFile (outfolder args ++ "/histograms.yoda") WriteMode $ \hand ->
    forM_ hists $ \(p, h) -> do
      let s = printYodaObj (T.pack p) . pure . H1DD $ over bins toArbBin h
      hPutStrLn hand $ T.unpack s
      hPutStrLn hand ""

  mapM_ tfileClose files

  where
    scaleWgt w e = e { eWeight = eWeight e * w }

    linesP fn = do
      s <- liftIO $ readFile fn
      each (lines s) >-> P.filter (not . null)

    readEvents isData f = do
      t <- ttree f "nominal_Loose"
      evalStateP t $ each [0..] >-> pipeTTree (readEvent isData) >-> P.concat

    readSumWeights f = do
      t <- ttree f "sumWeights"
      evalStateP t
        $ each [0..]
          >-> pipeTTree (float2Double <$> readBranch "totalEventsWeighted")

    readDSID f = do
      t <- ttree f "sumWeights"
      evalStateP t $ each [0..] >-> pipeTTree (readBranch "dsid")
