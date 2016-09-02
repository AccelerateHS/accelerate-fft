{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.CUDA
-- Copyright   : [2016] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
--

module Data.Array.Accelerate.Math.FFT.CUDA (

  fft1D,
  fft2D,
  fft3D,

) where

import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Twine
import Data.Array.Accelerate.Data.Complex

import Data.Array.Accelerate.CUDA.Foreign
import Data.Array.Accelerate.Array.Sugar                            as S hiding ( allocateArray )
import Data.Array.Accelerate.Type

import Foreign.Storable
import Foreign.CUDA.Analysis
import qualified Foreign.CUDA.FFT                                   as FFT
import qualified Foreign.CUDA.Driver                                as CUDA hiding ( device )
import qualified Foreign.CUDA.Driver.Context                        as CUDA ( device )

import Control.Concurrent.MVar
import Control.Exception
import Control.Monad
import Data.Maybe
import System.IO.Unsafe


fft1D :: IsFloating e
      => Mode
      -> CUDAForeignAcc (Vector (Complex e) -> (Vector (Complex e)))
fft1D mode = CUDAForeignAcc "fft1D" $ liftAtoC (cuFFT mode)

fft2D :: IsFloating e
      => Mode
      -> CUDAForeignAcc (Array DIM2 (Complex e) -> (Array DIM2 (Complex e)))
fft2D mode = CUDAForeignAcc "fft2D" $ liftAtoC (cuFFT mode)

fft3D :: IsFloating e
      => Mode
      -> CUDAForeignAcc (Array DIM3 (Complex e) -> (Array DIM3 (Complex e)))
fft3D mode = CUDAForeignAcc "fft3D" $ liftAtoC (cuFFT mode)


liftAtoC
    :: forall sh e. (Shape sh, IsFloating e)
    => (Stream -> Array (sh:.Int) e -> CIO (Array (sh:.Int) e))
    -> Stream
    -> Array (sh:.Int) (Complex e)
    -> CIO (Array (sh:.Int) (Complex e))
liftAtoC f s =
  case floatingType :: FloatingType e of
    TypeFloat{}   -> c2a s <=< f s <=< a2c s
    TypeDouble{}  -> c2a s <=< f s <=< a2c s
    TypeCFloat{}  -> c2a s <=< f s <=< a2c s
    TypeCDouble{} -> c2a s <=< f s <=< a2c s


-- | Call the cuFFT library to execute the FFT (inplace)
--
cuFFT :: forall sh e. (Shape sh, IsFloating e)
      => Mode
      -> Stream
      -> Array (sh:.Int) e
      -> CIO (Array (sh:.Int) e)
cuFFT mode st arr =
  withScalarArrayPtr arr st $ \d_arr -> liftIO $ do
    let sh :. sz = shape arr
    p <- plan (sh :. sz `quot` 2) (undefined::e)  -- recall this is an array of packed (Vec2 e)
    FFT.setStream p st
    case floatingType :: FloatingType e of
      TypeFloat{}   -> FFT.execC2C p d_arr d_arr (signOfMode mode) >> return arr
      TypeDouble{}  -> FFT.execZ2Z p d_arr d_arr (signOfMode mode) >> return arr
      TypeCFloat{}  -> FFT.execC2C p d_arr d_arr (signOfMode mode) >> return arr
      TypeCDouble{} -> FFT.execZ2Z p d_arr d_arr (signOfMode mode) >> return arr


-- | Convert an unzipped Accelerate array of complex numbers into a (new) packed
-- array suitable for use with CUFFT.
--
a2c :: forall sh e. (Shape sh, Elt e, IsFloating e, Storable (DevicePtrs e))
    => Stream
    -> Array (sh:.Int) (Complex e)
    -> CIO (Array (sh:.Int) e)              -- this is really a packed array of (Vec2 e) type
a2c st arr | FloatingDict <- floatingDict (floatingType :: FloatingType e) = do
  let
      sh :. sz  = shape arr
      n         = size sh * sz
  --
  cs <- allocateArray (sh :. 2*sz)
  withComplexArrayPtrs arr st $ \d_re d_im -> do
  withScalarArrayPtr   cs  st $ \d_cs      -> liftIO $ do
    mdl  <- twine (sizeOf (undefined::e))
    pack <- CUDA.getFun mdl "interleave"
    dev  <- CUDA.device
    prp  <- CUDA.props dev
    regs <- CUDA.requires pack CUDA.NumRegs
    let
        blockSize = 256
        sharedMem = 0
        maxBlocks = maxResidentBlocks prp blockSize regs sharedMem
        numBlocks = maxBlocks `min` ((n + blockSize - 1) `div` blockSize)
    --
    CUDA.launchKernel pack (numBlocks,1,1) (blockSize,1,1) sharedMem (Just st)
      [ CUDA.VArg d_cs, CUDA.VArg d_re, CUDA.VArg d_im, CUDA.IArg (fromIntegral n) ]
    return cs


-- | Convert a packed array of complex numbers into a (new) unzipped Accelerate
-- array.
--
c2a :: forall sh e. (Shape sh, Elt e, IsFloating e, Storable (DevicePtrs e))
    => Stream
    -> Array (sh:.Int) e
    -> CIO (Array (sh:.Int) (Complex e))
c2a st cs | FloatingDict <- floatingDict (floatingType :: FloatingType e) = do
  let
      sh :. sz2 = shape cs
      sz        = sz2 `quot` 2
      n         = size sh * sz
  --
  arr <- allocateArray (sh :. sz)
  withComplexArrayPtrs arr st $ \d_re d_im -> do
  withScalarArrayPtr   cs  st $ \d_cs      -> liftIO $ do
    mdl    <- twine (sizeOf (undefined::e))
    unpack <- CUDA.getFun mdl "deinterleave"
    dev    <- CUDA.device
    prp    <- CUDA.props dev
    regs   <- CUDA.requires unpack CUDA.NumRegs
    let
        blockSize = 256
        sharedMem = 0
        maxBlocks = maxResidentBlocks prp blockSize regs sharedMem
        numBlocks = maxBlocks `min` ((n + blockSize - 1) `div` blockSize)
    --
    CUDA.launchKernel unpack (numBlocks,1,1) (blockSize,1,1) sharedMem (Just st)
      [ CUDA.VArg d_re, CUDA.VArg d_im, CUDA.VArg d_cs, CUDA.IArg (fromIntegral n) ]
    return arr


-- | Generate an execute plan for a given type and size of FFT. These plans are
-- cached so that subsequent invocations are quicker.
--
plan :: forall sh e. (Shape sh, IsFloating e) => sh -> e -> IO FFT.Handle
plan (shapeToList -> sh) _ =
  modifyMVar fft_plans $ \ps ->
    case lookup (ty, sh) ps of
      Just p  -> return (ps, p)
      Nothing -> do
        p <- case sh of
               [w]     -> FFT.plan1D     w ty 1
               [w,h]   -> FFT.plan2D   h w ty
               [w,h,d] -> FFT.plan3D d h w ty
               _       -> error "cuFFT only supports 1D, 2D, and 3D transforms"
        return (((ty,sh),p) : ps, p)
  where
    ty = case floatingType :: FloatingType e of
           TypeFloat{}   -> FFT.C2C
           TypeDouble{}  -> FFT.Z2Z
           TypeCFloat{}  -> FFT.C2C
           TypeCDouble{} -> FFT.Z2Z


-- | Load the module to convert between SoA and AoS representation for the given
-- type. This is cached for subsequent reuse.
--
twine :: Int -> IO CUDA.Module
twine bitsize = do
  ctx <- fromMaybe (error "could not determine current CUDA context") `fmap` CUDA.get
  modifyMVar ptx_twine_modules $ \ms -> do
    case lookup (bitsize,ctx) ms of
      Just m  -> return (ms, m)
      Nothing -> do
        m <- CUDA.loadData $ case bitsize of
                               4 -> ptx_twine_f32
                               8 -> ptx_twine_f64
                               _ -> error "cuFFT only supports Float and Double"
        return (((bitsize,ctx), m) : ms, m)


-- | Dig out the two device pointers for an unzipped array of complex numbers.
--
withComplexArrayPtrs
    :: forall sh e a. IsFloating e
    => Array sh (Complex e)
    -> Stream
    -> (DevicePtrs e -> DevicePtrs e -> CIO a)
    -> CIO a
withComplexArrayPtrs arr st k
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withDevicePtrs arr (Just st) $ \(((),p1),p2) -> k p1 p2
      TypeDouble{}  -> withDevicePtrs arr (Just st) $ \(((),p1),p2) -> k p1 p2
      TypeCDouble{} -> withDevicePtrs arr (Just st) $ \(((),p1),p2) -> k p1 p2
      TypeCFloat{}  -> withDevicePtrs arr (Just st) $ \(((),p1),p2) -> k p1 p2

-- | Dig out the device pointer for a scalar array
--
withScalarArrayPtr
    :: forall sh e a. IsFloating e
    => Array sh e
    -> Stream
    -> (DevicePtrs e -> CIO a)
    -> CIO a
withScalarArrayPtr arr st k
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withDevicePtrs arr (Just st) $ \p -> k p
      TypeDouble{}  -> withDevicePtrs arr (Just st) $ \p -> k p
      TypeCDouble{} -> withDevicePtrs arr (Just st) $ \p -> k p
      TypeCFloat{}  -> withDevicePtrs arr (Just st) $ \p -> k p


-- Cache the FFT planning step for faster repeat evaluations.
{-# NOINLINE fft_plans #-}
fft_plans :: MVar [((FFT.Type, [Int]), FFT.Handle)]
fft_plans = unsafePerformIO $ do
  mv <- newMVar []
  _  <- mkWeakMVar mv
      $ withMVar mv
      $ mapM_ (\(_,p) -> FFT.destroy p)
  return mv

-- Cache the functions which convert between SoA and AoS format.
{-# NOINLINE ptx_twine_modules #-}
ptx_twine_modules :: MVar [((Int, CUDA.Context), CUDA.Module)]
ptx_twine_modules = unsafePerformIO $ do
  mv <- newMVar []
  _  <- mkWeakMVar mv
      $ withMVar mv
      $ mapM_ (\((_,ctx),mdl) -> bracket_ (CUDA.push ctx) CUDA.pop (CUDA.unload mdl))
  return mv

