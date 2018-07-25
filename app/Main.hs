{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main where

import           Atlas
import           Atlas.CrossSections
import           Control.Applicative        (ZipList (..))
import           Control.Comonad            (duplicate, extract)
import qualified Control.Foldl              as F
import           Control.Lens               (view)
import           Control.Monad              (forM_)
import           Control.Monad.Writer.Class
import           Data.Functor.Compose
import           Data.Maybe                 (isJust)
import           Data.Monoid                ((<>))
import           Data.TFile
import           Data.TTree
import           Debug.Trace
import           GHC.Float
import           Options.Generic
import           Pipes
import           Pipes.Lift
import qualified Pipes.Prelude              as P
import           System.IO

data Args =
  Args
  { outfile  :: String
  , infiles  :: String
  , xsecfile :: String
  , nevents  :: Maybe Int
  } deriving (Show, Generic)

instance ParseRecord Args where


data LargeJet =
  LargeJet
  { ljFourMom :: PtEtaPhiE
  , ljMass    :: Double
  } deriving Show


readLargeJets :: (MonadIO m, MonadThrow m) => TreeRead m [LargeJet]
readLargeJets = do
  pts <- fmap ((/1e3) . float2Double) <$> readBranch "ljet_pt"
  etas <- fmap float2Double <$> readBranch "ljet_eta"
  phis <- fmap float2Double <$> readBranch "ljet_phi"
  es <- fmap ((/1e3) . float2Double) <$> readBranch "ljet_e"
  ms <- fmap ((/1e3) . float2Double) <$> readBranch "ljet_m"

  let ljFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ LargeJet <$> ljFourMoms <*> ms


newtype RCJet =
  RCJet
  { rcFourMom :: PtEtaPhiE
  } deriving Show


readRCJets :: (MonadIO m, MonadThrow m) => TreeRead m [RCJet]
readRCJets = do
  pts <- fmap ((/1e3) . float2Double) <$> readBranch "rcjet_pt"
  etas <- fmap float2Double <$> readBranch "rcjet_eta"
  phis <- fmap float2Double <$> readBranch "rcjet_phi"
  es <- fmap ((/1e3) . float2Double) <$> readBranch "rcjet_e"

  let rcFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ RCJet <$> rcFourMoms


data Jet =
  Jet
  { jFourMom :: PtEtaPhiE
  , jDL177   :: Bool
  } deriving Show


readJets :: (MonadIO m, MonadThrow m) => TreeRead m [Jet]
readJets = do
  pts <- fmap ((/1e3) . float2Double) <$> readBranch "jet_pt"
  etas <- fmap float2Double <$> readBranch "jet_eta"
  phis <- fmap float2Double <$> readBranch "jet_phi"
  es <- fmap ((/1e3) . float2Double) <$> readBranch "jet_e"
  btags <- fmap cToB <$> readBranch "jet_isbtagged_DL1_77"

  let ljFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ Jet <$> ljFourMoms <*> btags

  where
    cToB :: CChar -> Bool
    cToB 0 = False
    cToB _ = True


readEvent :: (MonadIO m, MonadThrow m) => TreeRead m (Maybe Double)
readEvent = do
  -- wgt <- float2Double <$> readBranch "weight_mc"
  jets <- readJets
  topJets <- take 2 . filter topTagged <$> readLargeJets

  let bjets = filter bTagged jets
      tjp4s = ljFourMom <$> topJets
      jets' = removeOverlap tjp4s jets

  if length tjp4s < 2 || length bjets /= 2 || length jets' < 3
    then return Nothing
    else
      let (tj1:tj2:_) = tjp4s
      in return . Just . view lvM $ tj1 <> tj2

  where
    removeOverlap ljp4s =
      filter (\j -> all (\lj -> lvDRRap lj (jFourMom j) > 1.0) ljp4s)

    bTagged (Jet _ tagged) = tagged

    topTagged (LargeJet p4 m) = view lvPt p4 > 200 && m > 80


main :: IO ()
main = do
  args <- getRecord "run-4t" :: IO Args
  xsec <- readXSecFile (xsecfile args)

  hSetBuffering stdout LineBuffering

  let filesP = linesP $ infiles args

  -- sow <-
  --   F.purely P.fold F.sum
  --   $ for filesP readSumWeights

  -- putStrLn $ "sum of weights: " ++ show sow

  withFile (outfile args) WriteMode $ \h -> do
    hSetBuffering h LineBuffering
    runEffect $ for filesP readTree >-> P.toHandle h

  where
    linesP fn = do
      s <- liftIO $ readFile fn
      each (lines s) >-> P.filter (not . null)

    readTree fn = do
      liftIO . putStrLn $ "running over file " ++ fn
      f <- tfileOpen fn
      t <- ttree f "nominal_Loose"
      evalStateP t
        $ each [0..]
          >-> pipeTTree readEvent
          >-> P.concat
          >-> P.map (\m -> "1.0, " ++ show m)

      tfileClose f

    readSumWeights fn = do
      liftIO . putStrLn $ "checking weights of file " ++ fn
      f <- tfileOpen fn
      t <- ttree f "sumWeights"
      evalStateP t
        $ each [0..]
          >-> pipeTTree (float2Double <$> readBranch "totalEventsWeighted")
