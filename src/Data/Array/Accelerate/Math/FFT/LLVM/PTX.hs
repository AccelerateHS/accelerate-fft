{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
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
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Lifetime

import Data.Array.Accelerate.LLVM.PTX.Foreign

import Foreign.CUDA.Ptr                                             ( DevicePtr, castDevPtr )
import qualified Foreign.CUDA.FFT                                   as FFT

import Control.Monad.Reader
import Data.Hashable
import Data.Proxy
import System.IO.Unsafe


fft :: forall sh e. (Shape sh, Numeric e)
    => Mode
    -> ForeignAcc (Array (sh:.Int) (Complex e) -> Array (sh:.Int) (Complex e))
fft mode
  | Just Refl <- matchShapeType @sh @DIM0 = fft1D mode
  | Just Refl <- matchShapeType @sh @DIM1 = ForeignAcc "cuda.fft2.many" $ fft' fft2DMany_plans mode
  | Just Refl <- matchShapeType @sh @DIM2 = ForeignAcc "cuda.fft3.many" $ fft' fft3DMany_plans mode
  | otherwise = $internalError "fft" "only for 1D..3D inner-dimension transforms"

fft1D :: Numeric e
      => Mode
      -> ForeignAcc (Vector (Complex e) -> Vector (Complex e))
fft1D mode = ForeignAcc "cuda.fft1d" $ fft' fft1D_plans mode

fft2D :: Numeric e
      => Mode
      -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
fft2D mode = ForeignAcc "cuda.fft2d" $ fft' fft2D_plans mode

fft3D :: Numeric e
      => Mode
      -> ForeignAcc (Array DIM3 (Complex e) -> Array DIM3 (Complex e))
fft3D mode = ForeignAcc "cuda.fft3d" $ fft' fft3D_plans mode


-- Internals
-- ---------

{-# INLINEABLE fft' #-}
fft' :: forall sh e. (Shape sh, Numeric e)
     => Plans (sh, FFT.Type)
     -> Mode
     -> Array sh (Complex e)
     -> Par PTX (Future (Array sh (Complex e)))
fft' plans mode =
  let
      go :: Numeric e => Array sh (Complex e) -> Par PTX (Future (Array sh (Complex e)))
      go ain = do
        let
            sh = shape ain
            t  = fftType (Proxy::Proxy e)
        --
        aout    <- allocateRemote sh
        stream  <- asks ptxStream
        future  <- new
        liftPar $
          withArray ain stream    $ \d_in  -> do
           withArray aout stream  $ \d_out -> do
            withPlan plans (sh,t) $ \h     -> do
              liftIO $ cuFFT (Proxy::Proxy e) h mode stream (castDevPtr d_in) (castDevPtr d_out)
        --
        put future aout
        return future
  in
  case numericR::NumericR e of
    NumericRfloat32 -> go
    NumericRfloat64 -> go


-- Execute the FFT
--
{-# INLINE cuFFT #-}
cuFFT :: forall e. Numeric e
      => Proxy e
      -> FFT.Handle
      -> Mode
      -> Stream
      -> DevicePtr (Complex e)
      -> DevicePtr (Complex e)
      -> IO ()
cuFFT _ p mode stream d_in d_out =
  withLifetime stream $ \s -> do
    FFT.setStream p s
    case numericR::NumericR e of
      NumericRfloat32 -> FFT.execC2C p (fftMode mode) d_in d_out
      NumericRfloat64 -> FFT.execZ2Z p (fftMode mode) d_in d_out

fftType :: forall e. Numeric e => Proxy e -> FFT.Type
fftType _ =
  case numericR::NumericR e of
    NumericRfloat32 -> FFT.C2C
    NumericRfloat64 -> FFT.Z2Z

fftMode :: Mode -> FFT.Mode
fftMode Forward = FFT.Forward
fftMode _       = FFT.Inverse


-- Plan caches
-- -----------

{-# NOINLINE fft1D_plans #-}
fft1D_plans :: Plans (DIM1, FFT.Type)
fft1D_plans
  = unsafePerformIO
  $ createPlan (\(Z:.n, t) -> FFT.plan1D n t 1)
               (\(Z:.n, t) -> fromEnum t `hashWithSalt` n)

{-# NOINLINE fft2D_plans #-}
fft2D_plans :: Plans (DIM2, FFT.Type)
fft2D_plans
  = unsafePerformIO
  $ createPlan (\(Z:.h:.w, t) -> FFT.plan2D h w t)
               (\(Z:.h:.w, t) -> fromEnum t `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft3D_plans #-}
fft3D_plans :: Plans (DIM3, FFT.Type)
fft3D_plans
  = unsafePerformIO
  $ createPlan (\(Z:.d:.h:.w, t) -> FFT.plan3D d h w t)
               (\(Z:.d:.h:.w, t) -> fromEnum t `hashWithSalt` d `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft2DMany_plans #-}
fft2DMany_plans :: Plans (DIM2, FFT.Type)
fft2DMany_plans
  = unsafePerformIO
  $ createPlan (\(Z:.h:.w, t) -> FFT.planMany [w] Nothing Nothing t h)
               (\(Z:.h:.w, t) -> fromEnum t `hashWithSalt` h `hashWithSalt` w)

{-# NOINLINE fft3DMany_plans #-}
fft3DMany_plans :: Plans (DIM3, FFT.Type)
fft3DMany_plans
  = unsafePerformIO
  $ createPlan (\(Z:.d:.h:.w, t) -> FFT.planMany [w] Nothing Nothing t (d*h))
               (\(Z:.d:.h:.w, t) -> fromEnum t `hashWithSalt` d `hashWithSalt` h `hashWithSalt` w)

