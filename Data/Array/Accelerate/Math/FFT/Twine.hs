{-# LANGUAGE ScopedTypeVariables #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Twine
-- Copyright   : [2016] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
--

module Data.Array.Accelerate.Math.FFT.Twine
  where

import Data.Array.Accelerate                                      as A
import Data.Array.Accelerate.Data.Complex


-- Interleave the real and imaginary components in a complex array and produce a
-- flattened vector. This allows us to mimic the array-of-struct representation
-- commonly used by FFT libraries to store complex numbers (CUFFT, FFTW).
--
-- We would really prefer to implement this with a zipWith of the two arrays,
-- but we can't represent the packed structure in Accelerate.
--
{-# NOINLINE interleave #-}
interleave :: Elt e => Acc (Vector (Complex e)) -> Acc (Vector e)
interleave arr = generate sh swizzle
  where
    reals       = A.map real arr
    imags       = A.map imag arr
    --
    sh          = index1 (2 * A.size arr)
    swizzle ix  =
      let i     = indexHead ix
          (j,k) = i `quotRem` 2
      in
      k ==* 0 ? ( reals A.!! j, imags A.!! j )


-- Deinterleave a vector into a complex array. Requires the array to have an
-- even number of elements.
--
{-# NOINLINE deinterleave #-}
deinterleave :: forall e. Elt e => Acc (Vector e) -> Acc (Vector (Complex e))
deinterleave arr = generate sh swizzle
  where
    sh         = index1 (A.size arr `quot` 2)
    swizzle ix =
      let i = indexHead ix `quot` 2
      in  lift ( arr A.!! i :+ arr A.!! (i+1) ) :: Exp (Complex e)


{-# RULES
  "interleave/deinterleave" forall x. deinterleave (interleave x) = x;
  "deinterleave/interleave" forall x. interleave (deinterleave x) = x
 #-}

