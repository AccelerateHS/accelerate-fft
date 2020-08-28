-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Mode
-- Copyright   : [2012..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
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

