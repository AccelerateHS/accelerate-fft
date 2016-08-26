-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Mode
-- Copyright   : [2012..2016] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
--               [2013..2014] Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.Mode
  where


data Mode = Forward | Reverse | Inverse
  deriving (Eq, Show)

signOfMode :: Num a => Mode -> a
signOfMode m
  = case m of
      Forward   -> -1
      Reverse   ->  1
      Inverse   ->  1

