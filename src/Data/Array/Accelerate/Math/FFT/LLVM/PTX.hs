{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections       #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX
-- Copyright   : [2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX (

  fft,
  fft1D,
  fft2D,
  fft3D,

) where

import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Type
import Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base
import Data.Array.Accelerate.Math.FFT.LLVM.PTX.Plans

import Data.Array.Accelerate.Analysis.Match
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.Representation.Array
import Data.Array.Accelerate.Representation.Shape
import Data.Array.Accelerate.Sugar.Elt
import Data.Primitive.Vec

import Data.Array.Accelerate.LLVM.PTX.Foreign

import Foreign.CUDA.Ptr                                             ( DevicePtr, castDevPtr )
import qualified Foreign.CUDA.FFT                                   as FFT

import Control.Monad.Reader
import Data.Hashable
import System.IO.Unsafe


fft :: forall sh e. HasCallStack
    => Mode
    -> ShapeR sh
    -> NumericR e
    -> ForeignAcc (Array sh (Vec2 e) -> Array sh (Vec2 e))
fft mode shR eR
  | Just Refl <- matchShapeR shR dim1 = fft1D mode eR
  | Just Refl <- matchShapeR shR dim2 = ForeignAcc "cuda.fft2.many" $ fft' fft2DMany_plans mode shR eR
  | Just Refl <- matchShapeR shR dim3 = ForeignAcc "cuda.fft3.many" $ fft' fft3DMany_plans mode shR eR
  | otherwise = internalError "only for 1D..3D inner-dimension transforms"

fft1D :: Mode -> NumericR e -> ForeignAcc (Array DIM1 (Vec2 e) -> Array DIM1 (Vec2 e))
fft1D mode eR = ForeignAcc "cuda.fft1d" $ fft' fft1D_plans mode dim1 eR

fft2D :: Mode -> NumericR e -> ForeignAcc (Array DIM2 (Vec2 e) -> Array DIM2 (Vec2 e))
fft2D mode eR = ForeignAcc "cuda.fft2d" $ fft' fft2D_plans mode dim2 eR

fft3D :: Mode -> NumericR e -> ForeignAcc (Array DIM3 (Vec2 e) -> Array DIM3 (Vec2 e))
fft3D mode eR = ForeignAcc "cuda.fft3d" $ fft' fft3D_plans mode dim3 eR


-- Internals
-- ---------

{-# INLINEABLE fft' #-}
fft' :: forall sh e.
        Plans (sh, FFT.Type)
     -> Mode
     -> ShapeR sh
     -> NumericR e
     -> Array sh (Vec2 e)
     -> Par PTX (Future (Array sh (Vec2 e)))
fft' plans mode shR eR =
  let
      go :: ArrayR (Array sh (Vec2 e)) -> Array sh (Vec2 e) -> Par PTX (Future (Array sh (Vec2 e)))
      go aR ain = do
        let
            sh = shape ain
            t  = fftType eR
        --
        aout    <- allocateRemote aR sh
        stream  <- asks ptxStream
        future  <- new
        liftPar $
          withArray eR ain stream   $ \d_in  -> do
           withArray eR aout stream $ \d_out -> do
            withPlan plans (sh,t)   $ \h     -> do
              liftIO $ cuFFT eR h mode stream (castDevPtr d_in) (castDevPtr d_out)
        --
        put future aout
        return future
  in
  case eR of
    NumericRfloat32 -> go (ArrayR shR (eltR @(Complex Float)))
    NumericRfloat64 -> go (ArrayR shR (eltR @(Complex Double)))

--
-- Execute the FFT
--
{-# INLINE cuFFT #-}
cuFFT :: NumericR e
      -> FFT.Handle
      -> Mode
      -> Stream
      -> DevicePtr (Complex e)
      -> DevicePtr (Complex e)
      -> IO ()
cuFFT eR p mode stream d_in d_out =
  withLifetime stream $ \s -> do
    FFT.setStream p s
    case eR of
      NumericRfloat32 -> FFT.execC2C p (fftMode mode) d_in d_out
      NumericRfloat64 -> FFT.execZ2Z p (fftMode mode) d_in d_out

fftType :: NumericR e -> FFT.Type
fftType NumericRfloat32 = FFT.C2C
fftType NumericRfloat64 = FFT.Z2Z

fftMode :: Mode -> FFT.Mode
fftMode Forward = FFT.Forward
fftMode _       = FFT.Inverse

-- Plan caches
-- -----------

{-# NOINLINE fft1D_plans #-}
fft1D_plans :: Plans (DIM1, FFT.Type)
fft1D_plans
  = unsafePerformIO
  $ createPlan (\(((), n), t) -> FFT.plan1D n t 1)
               (\(((), n), t) -> fromEnum t `hashWithSalt` n)

{-# NOINLINE fft2D_plans #-}
fft2D_plans :: Plans (DIM2, FFT.Type)
fft2D_plans
  = unsafePerformIO
  $ createPlan (\((((),h),w), t) -> FFT.plan2D h w t)
               (\((((),h),w), t) -> fromEnum t `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft3D_plans #-}
fft3D_plans :: Plans (DIM3, FFT.Type)
fft3D_plans
  = unsafePerformIO
  $ createPlan (\(((((),d),h),w), t) -> FFT.plan3D d h w t)
               (\(((((),d),h),w), t) -> fromEnum t `hashWithSalt` d `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft2DMany_plans #-}
fft2DMany_plans :: Plans (DIM2, FFT.Type)
fft2DMany_plans
  = unsafePerformIO
  $ createPlan (\((((),h),w), t) -> FFT.planMany [w] Nothing Nothing t h)
               (\((((),h),w), t) -> fromEnum t `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft3DMany_plans #-}
fft3DMany_plans :: Plans (DIM3, FFT.Type)
fft3DMany_plans
  = unsafePerformIO
  $ createPlan (\(((((()),d),h),w), t) -> FFT.planMany [w] Nothing Nothing t (d*h))
               (\(((((()),d),h),w), t) -> fromEnum t `hashWithSalt` d `hashWithSalt` h `hashWithSalt` w)

