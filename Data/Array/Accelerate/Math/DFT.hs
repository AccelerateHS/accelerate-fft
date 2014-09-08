{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT
-- Copyright   : [2012] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Manuel M T Chakravarty <chak@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Compute the Discrete Fourier Transform (DFT) along the lower order dimension
-- of an array.
--
-- This uses a naïve algorithm which takes O(n^2) time. However, you can
-- transform an array with an arbitrary extent, unlike with FFT which requires
-- each dimension to be a power of two.
--
-- The `dft` and `idft` functions compute the roots of unity as needed. If you
-- need to transform several arrays with the same extent than it is faster to
-- compute the roots once using `rootsOfUnity` or `inverseRootsOfUnity`
-- respectively, then call `dftG` directly.
--
-- You can also compute single values of the transform using `dftGS`
--
module Data.Array.Accelerate.Math.DFT (

  dft, idft, dftG, dftGS,

) where

import Prelude                                  as P hiding ((!!))
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Math.DFT.Roots
import Data.Array.Accelerate.Data.Complex


-- | Compute the DFT along the low order dimension of an array
--
dft :: (Shape sh, Slice sh, Elt e, IsFloating e)
    => Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
dft v = dftG (rootsOfUnity (shape v)) v


-- | Compute the inverse DFT along the low order dimension of an array
--
idft :: (Shape sh, Slice sh, Elt e, IsFloating e)
     => Acc (Array (sh:.Int) (Complex e))
     -> Acc (Array (sh:.Int) (Complex e))
idft v
  = let sh      = shape v
        n       = indexHead sh
        roots   = inverseRootsOfUnity sh
        scale   = lift (A.fromIntegral n :+ constant 0)
    in
    A.map (/scale) $ dftG roots v


-- | Generic function for computation of forward and inverse DFT. This function
--   is also useful if you transform many arrays of the same extent, and don't
--   want to recompute the roots for each one.
--
--   The extent of the input and roots must match.
--
dftG :: forall sh e. (Shape sh, Slice sh, Elt e, IsFloating e)
     => Acc (Array (sh:.Int) (Complex e))       -- ^ roots of unity
     -> Acc (Array (sh:.Int) (Complex e))       -- ^ input array
     -> Acc (Array (sh:.Int) (Complex e))
dftG roots arr
  = A.fold (+) (constant (0 :+ 0))
  $ A.zipWith (*) arr' roots'
  where
    base        = shape arr
    l           = indexHead base
    extend      = lift (base :. shapeSize base)

    -- Extend the entirety of the input arrays into a higher dimension, reading
    -- roots from the appropriate places and then reduce along this axis.
    --
    -- In the calculation for 'roots'', 'i' is the index into the extended
    -- dimension, with corresponding base index 'ix' which we are attempting to
    -- calculate the single DFT value of. The rest proceeds as per 'dftGS'.
    --
    arr'        = A.generate extend (\ix' -> let i = indexHead ix' in arr !! i)
    roots'      = A.generate extend (\ix' -> let ix :. i    = unlift ix'
                                                 sh :. n    = unlift (fromIndex base i) :: Exp sh :. Exp Int
                                                 k          = indexHead ix
                                             in
                                             roots ! lift (sh :. (k*n) `mod` l))


-- | Compute a single value of the DFT.
--
dftGS :: forall sh e. (Shape sh, Slice sh, Elt e, IsFloating e)
      => Exp (sh :. Int)                        -- ^ index of the value we want
      -> Acc (Array (sh:.Int) (Complex e))      -- ^ roots of unity
      -> Acc (Array (sh:.Int) (Complex e))      -- ^ input array
      -> Acc (Scalar (Complex e))
dftGS ix roots arr
  = let k = indexHead ix
        l = indexHead (shape arr)

        -- all the roots we need to multiply with
        roots'  = A.generate (shape arr)
                             (\ix' -> let sh :. n = unlift ix'  :: Exp sh :. Exp Int
                                      in  roots ! lift (sh :. (k*n) `mod` l))
    in
    A.foldAll (+) (constant (0 :+ 0)) $ A.zipWith (*) arr roots'

