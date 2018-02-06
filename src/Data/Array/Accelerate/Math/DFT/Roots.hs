{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT.Roots
-- Copyright   : [2012..2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
--               [2013..2017] Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
module Data.Array.Accelerate.Math.DFT.Roots (

  rootsOfUnity, inverseRootsOfUnity,

) where

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Data.Complex


-- | Calculate the roots of unity for the forward transform
--
rootsOfUnity
    :: (Shape sh, Slice sh, Elt (Complex e), A.Floating e, A.FromIntegral Int e)
    => Exp (sh :. Int)
    -> Acc (Array (sh:.Int) (Complex e))
rootsOfUnity sh =
  let n = A.fromIntegral (A.indexHead sh)
  in
  A.generate sh (\ix -> let i = A.fromIntegral (A.indexHead ix)
                            k = 2 * pi * i / n
                        in
                        A.lift ( cos k :+ (-sin k) ))


-- | Calculate the roots of unity for an inverse transform
--
inverseRootsOfUnity
    :: (Shape sh, Slice sh, Elt (Complex e), A.Floating e, A.FromIntegral Int e)
    => Exp (sh :. Int)
    -> Acc (Array (sh:.Int) (Complex e))
inverseRootsOfUnity sh =
  let n = A.fromIntegral (A.indexHead sh)
  in
  A.generate sh (\ix -> let i = A.fromIntegral (A.indexHead ix)
                            k = 2 * pi * i / n
                        in
                        A.lift ( cos k :+ sin k ))

