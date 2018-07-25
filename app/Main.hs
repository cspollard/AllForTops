{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}

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
-- import           Data.Histogram.Bin.BinF
import           Data.Monoid                ((<>))
import           Data.TFile
import           Data.TTree
import           GHC.Float
import           Options.Generic
import           Pipes
import           Pipes.Lift
import qualified Pipes.Prelude              as P

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
  pts <- fmap float2Double <$> readBranch "ljet_pt"
  etas <- fmap float2Double <$> readBranch "ljet_eta"
  phis <- fmap float2Double <$> readBranch "ljet_phi"
  es <- fmap float2Double <$> readBranch "ljet_e"
  ms <- fmap float2Double <$> readBranch "ljet_m"

  let ljFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ LargeJet <$> ljFourMoms <*> ms


newtype RCJet =
  RCJet
  { rcFourMom :: PtEtaPhiE
  } deriving Show

readRCJets :: (MonadIO m, MonadThrow m) => TreeRead m [RCJet]
readRCJets = do
  pts <- fmap float2Double <$> readBranch "rcjet_pt"
  etas <- fmap float2Double <$> readBranch "rcjet_eta"
  phis <- fmap float2Double <$> readBranch "rcjet_phi"
  es <- fmap float2Double <$> readBranch "rcjet_e"

  let rcFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ RCJet <$> rcFourMoms


data Jet =
  Jet
  { jFourMom :: PtEtaPhiE
  , jDL177   :: Bool
  } deriving Show


readJets :: (MonadIO m, MonadThrow m) => TreeRead m [Jet]
readJets = do
  pts <- fmap float2Double <$> readBranch "jet_pt"
  etas <- fmap float2Double <$> readBranch "jet_eta"
  phis <- fmap float2Double <$> readBranch "jet_phi"
  es <- fmap float2Double <$> readBranch "jet_e"
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
  js <- readJets
  rcjs <- readRCJets
  let nbjs = foldr (\(Jet _ tagged) s -> if tagged then s+1 else s) 0 js
      js' = removeOverlap rcjs js
      rcp4s = take 2 . filter ((> 80000) . view lvM) $ rcFourMom <$> rcjs

  if nbjs < 3 || length js' < 4 || length rcjs < 2
    then return Nothing
    else
      let (rc1:rc2:_) = rcp4s
      in return . Just . view lvM $ rc1 <> rc2

  where
    removeOverlap rcjs js =
      let rcp4s = rcFourMom <$> rcjs
      in filter (\j -> all (\rc -> lvDRRap rc (jFourMom j) > 1.0) rcp4s) js


main :: IO ()
main = do
  args <- getRecord "run-4t" :: IO Args
  xsec <- readXSecFile (xsecfile args)


  runEffect
    $ for (linesP $ infiles args) readTree >-> P.mapM_ printLine


  where
    linesP fn = do
      s <- liftIO $ readFile fn
      each (lines s) >-> P.filter (not . null)

    readTree fn = do
      f <- tfileOpen fn
      t <- ttree f "nominal_Loose"
      evalStateP t
        $ each [0..]
          >-> pipeTTree readEvent
          >-> P.concat
          >-> P.map (1.0,)

      tfileClose f

    printLine (w, m) = putStr $ show w ++ ", " ++ show m
