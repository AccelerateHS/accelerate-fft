-- |
-- Module      : TestPTX
-- Copyright   : [2017..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module TestPTX where

import Test.FFT
import Test.Tasty
import Data.Array.Accelerate.LLVM.PTX                            as PTX

main :: IO ()
main = defaultMain (testFFT PTX.runN)

