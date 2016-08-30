{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
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

import Data.FileEmbed
import Data.ByteString                                            ( ByteString )


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


-- Embedded PTX code for interleave and deinterleave for 32- and 64-bit floating
-- point numbers respectively. These can be loaded and executed by the CUDA
-- driver at runtime as required.
--
-- The PTX code was compiled for SM-2.0 and 64-bit address space (the default
-- settings of nvcc-7.5), but the code is simple enough that the CUDA device
-- driver should be able to compile it for the actual target architecture
-- without issue. This has been confirmed with respect to SM, but I don't have
-- a 32-bit machine available to test that aspect with.
--

ptx_twine_f32 :: ByteString
ptx_twine_f32 = $(makeRelativeToProject "cubits/twine_f32.ptx" >>= embedFile)

ptx_twine_f64 :: ByteString
ptx_twine_f64 = $(makeRelativeToProject "cubits/twine_f64.ptx" >>= embedFile)

