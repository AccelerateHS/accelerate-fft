{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.Native
-- Copyright   : [2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.Native (

  fft,
  fft1D,
  fft2D,
  fft3D,

) where

import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Type
import Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix
import Data.Array.Accelerate.Math.FFT.LLVM.Native.Base

import Data.Array.Accelerate
import Data.Array.Accelerate.Analysis.Match
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Error

import Data.Array.Accelerate.LLVM.Native.Foreign

import Data.Array.CArray                                            ( CArray )
import Math.FFT.Base                                                ( FFTWReal )
import Prelude                                                      as P
import qualified Math.FFT                                           as FFT


fft :: forall sh e. (Shape sh, Numeric e)
    => Mode
    -> ForeignAcc (Array sh (Complex e) -> Array sh (Complex e))
fft mode
  = ForeignAcc (nameOf mode (undefined::sh))
  $ case numericR::NumericR e of
      NumericRfloat32 -> go
      NumericRfloat64 -> go
  where
    go :: FFTWReal e => Array sh (Complex e) -> LLVM Native (Array sh (Complex e))
    go | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM1) = liftCtoA (FFT.dftGU (signOf mode) flags [0] `ix` (undefined :: (Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM2) = liftCtoA (FFT.dftGU (signOf mode) flags [1] `ix` (undefined :: (Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM3) = liftCtoA (FFT.dftGU (signOf mode) flags [2] `ix` (undefined :: (Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM4) = liftCtoA (FFT.dftGU (signOf mode) flags [3] `ix` (undefined :: (Int,Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM5) = liftCtoA (FFT.dftGU (signOf mode) flags [4] `ix` (undefined :: (Int,Int,Int,Int,Int)))
       | otherwise = $internalError "fft" "only for 1D..5D inner-dimension transforms"
    --
    ix :: (a i r -> a i r) -> i -> (a i r -> a i r)
    ix f _ = f


fft1D :: forall e. Numeric e
      => Mode
      -> ForeignAcc (Array DIM1 (Complex e) -> Array DIM1 (Complex e))
fft1D mode
  = ForeignAcc (nameOf mode (undefined::DIM1))
  $ case numericR::NumericR e of
      NumericRfloat32 -> liftCtoA go
      NumericRfloat64 -> liftCtoA go
  where
    go :: FFTWReal r => CArray Int (Complex r) -> CArray Int (Complex r)
    go = FFT.dftGU (signOf mode) flags [0]

fft2D :: forall e. Numeric e
      => Mode
      -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
fft2D mode
  = ForeignAcc (nameOf mode (undefined::DIM2))
  $ case numericR::NumericR e of
      NumericRfloat32 -> liftCtoA go
      NumericRfloat64 -> liftCtoA go
  where
    go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1]

fft3D :: forall e. Numeric e
      => Mode
      -> ForeignAcc (Array DIM3 (Complex e) -> Array DIM3 (Complex e))
fft3D mode
  = ForeignAcc (nameOf mode (undefined::DIM3))
  $ case numericR::NumericR e of
      NumericRfloat32 -> liftCtoA go
      NumericRfloat64 -> liftCtoA go
  where
    go :: FFTWReal r => CArray (Int,Int,Int) (Complex r) -> CArray (Int,Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1,2]


liftCtoA
    :: forall ix sh e. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix, Numeric e)
    => (CArray ix (Complex e) -> CArray ix (Complex e))
    -> Array sh (Complex e)
    -> LLVM Native (Array sh (Complex e))
liftCtoA f a =
  liftIO $ withCArray a (fromCArray . f)

