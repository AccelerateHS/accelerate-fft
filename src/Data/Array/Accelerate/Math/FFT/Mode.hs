-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Mode
-- Copyright   : [2012..2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
--               [2013..2017] Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.Mode
  where


data Mode
  = Forward         -- ^ Forward DFT
  | Reverse         -- ^ Inverse DFT, un-normalised
  | Inverse         -- ^ Inverse DFT, normalised
  deriving (Eq, Show)

signOfMode :: Num a => Mode -> a
signOfMode m
  = case m of
      Forward   -> -1
      Reverse   ->  1
      Inverse   ->  1

