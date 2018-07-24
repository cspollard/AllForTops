{-# LANGUAGE DeriveGeneric       #-}
{-# LANGUAGE OverloadedStrings   #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main where

import           Atlas
import           Atlas.CrossSections
import           Control.Applicative (ZipList (..))
import qualified Control.Foldl       as F
import           Control.Monad       (forM_)
import           Data.TFile
import           Data.TTree
import           GHC.Float
import           Options.Generic
import           Pipes
import           Pipes.Lift
import qualified Pipes.Prelude       as P

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


data RCJet =
  RCJet
  { rcFourMom :: PtEtaPhiE
  , rcMass    :: Double
  } deriving Show

readRCJets :: (MonadIO m, MonadThrow m) => TreeRead m [RCJet]
readRCJets = do
  pts <- fmap float2Double <$> readBranch "rcjet_pt"
  etas <- fmap float2Double <$> readBranch "rcjet_eta"
  phis <- fmap float2Double <$> readBranch "rcjet_phi"
  es <- fmap float2Double <$> readBranch "rcjet_e"
  ms <- fmap float2Double <$> readBranch "rcjet_m"

  let ljFourMoms = PtEtaPhiE <$> pts <*> etas <*> phis <*> es

  return . getZipList $ RCJet <$> ljFourMoms <*> ms


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


eventRead :: (MonadIO m, MonadThrow m) => TreeRead m ([Jet], [LargeJet])
eventRead = (,) <$> readJets <*> readLargeJets


main :: IO ()
main = do
  args <- getRecord "run-4t" :: IO Args
  xsec <- readXSecFile (xsecfile args)

  -- get the list of input trees
  fns <- filter (not . null) . lines <$> readFile (infiles args)

  forM_ fns
    $ \fn -> do
      f <- tfileOpen fn
      t <- ttree f "nominal_Loose"
      runEffect . evalStateP t
        $ each [0..]
          >-> pipeTTree eventRead
          >-> P.print

      tfileClose f
