{-# LANGUAGE TypeOperators #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT.Roots
-- Copyright   : [2012] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Manuel M T Chakravarty <chak@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
module Data.Array.Accelerate.Math.DFT.Roots (

  rofu, irofu,

) where

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Math.Complex


-- | Calculate the roots of unity for the forward transform
--
rofu :: (Elt e, IsFloating e, Shape sh, Slice sh)
     => Exp (sh :. Int)
     -> Acc (Array (sh:.Int) (Complex e))
rofu sh =
  let n = A.fromIntegral (A.indexHead sh)
  in
  A.generate sh (\ix -> let i = A.fromIntegral (A.indexHead ix)
                            k = 2 * pi * i / n
                        in
                        A.lift ( cos k, -sin k ))


-- | Calculate the roots of unity for an inverse transform
--
irofu :: (Elt e, IsFloating e, Shape sh, Slice sh)
      => Exp (sh :. Int)
      -> Acc (Array (sh:.Int) (Complex e))
irofu sh =
  let n = A.fromIntegral (A.indexHead sh)
  in
  A.generate sh (\ix -> let i = A.fromIntegral (A.indexHead ix)
                            k = 2 * pi * i / n
                        in
                        A.lift ( cos k, sin k ))

