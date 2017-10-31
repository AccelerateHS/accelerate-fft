-- |
-- Module      : TestNative
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module TestNative where

import Test.FFT
import Test.Tasty
import Data.Array.Accelerate.LLVM.Native                            as CPU

main :: IO ()
main = defaultMain (testFFT CPU.runN)

