{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
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

import Data.Array.Accelerate.Analysis.Match
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Representation.Array
import Data.Array.Accelerate.Representation.Shape
import Data.Array.Accelerate.Sugar.Elt

import Data.Primitive.Vec

import Data.Array.Accelerate.LLVM.Native.Foreign

import Data.Proxy
import Data.Array.CArray                                            ( CArray )
import Math.FFT.Base                                                ( FFTWReal )
import Prelude                                                      as P
import qualified Math.FFT                                           as FFT


fft :: forall sh e. HasCallStack
    => Mode
    -> ShapeR sh
    -> NumericR e
    -> ForeignAcc (Array sh (Vec2 e) -> Array sh (Vec2 e))
fft mode shR eR
  = ForeignAcc (nameOf mode shR)
  $ case eR of
      NumericRfloat32 -> go
      NumericRfloat64 -> go
  where
    go :: FFTWReal e => Array sh (Vec2 e) -> Par Native (Future (Array sh (Vec2 e)))
    go | Just Refl <- matchShapeR shR dim1 = liftCtoA shR eR (FFT.dftGU (signOf mode) flags [0] `ix` (Proxy :: Proxy (Int)))
       | Just Refl <- matchShapeR shR dim2 = liftCtoA shR eR (FFT.dftGU (signOf mode) flags [1] `ix` (Proxy :: Proxy (Int,Int)))
       | Just Refl <- matchShapeR shR dim3 = liftCtoA shR eR (FFT.dftGU (signOf mode) flags [2] `ix` (Proxy :: Proxy (Int,Int,Int)))
       | Just Refl <- matchShapeR shR dim4 = liftCtoA shR eR (FFT.dftGU (signOf mode) flags [3] `ix` (Proxy :: Proxy (Int,Int,Int,Int)))
       | Just Refl <- matchShapeR shR dim5 = liftCtoA shR eR (FFT.dftGU (signOf mode) flags [4] `ix` (Proxy :: Proxy (Int,Int,Int,Int,Int)))
       | otherwise = internalError "only for 1D..5D inner-dimension transforms"
    --
    ix :: (a i r -> a i r) -> proxy i -> (a i r -> a i r)
    ix f _ = f

    dim4 = ShapeRsnoc dim3
    dim5 = ShapeRsnoc dim4

fft1D :: Mode -> NumericR e -> ForeignAcc (Array DIM1 (Vec2 e) -> Array DIM1 (Vec2 e))
fft1D mode eR
  = ForeignAcc (nameOf mode dim1)
  $ case eR of
      NumericRfloat32 -> liftCtoA dim1 eR go
      NumericRfloat64 -> liftCtoA dim1 eR go
  where
    go :: FFTWReal r => CArray Int (Complex r) -> CArray Int (Complex r)
    go = FFT.dftGU (signOf mode) flags [0]

fft2D :: Mode -> NumericR e -> ForeignAcc (Array DIM2 (Vec2 e) -> Array DIM2 (Vec2 e))
fft2D mode eR
  = ForeignAcc (nameOf mode dim2)
  $ case eR of
      NumericRfloat32 -> liftCtoA dim2 eR go
      NumericRfloat64 -> liftCtoA dim2 eR go
  where
    go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1]

fft3D :: Mode -> NumericR e -> ForeignAcc (Array DIM3 (Vec2 e) -> Array DIM3 (Vec2 e))
fft3D mode eR
  = ForeignAcc (nameOf mode dim3)
  $ case eR of
      NumericRfloat32 -> liftCtoA dim3 eR go
      NumericRfloat64 -> liftCtoA dim3 eR go
  where
    go :: FFTWReal r => CArray (Int,Int,Int) (Complex r) -> CArray (Int,Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1,2]

liftCtoA
    :: forall ix sh e. (IxShapeR (EltR ix) ~ sh, Elt ix)
    => ShapeR sh
    -> NumericR e
    -> (CArray ix (Complex e) -> CArray ix (Complex e))
    -> Array sh (Vec2 e)
    -> Par Native (Future (Array sh (Vec2 e)))
liftCtoA shR eR f a =
  newFull =<< liftIO (withCArray shR eR a (fromCArray shR eR . f))

