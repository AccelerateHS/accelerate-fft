{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}
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

module Data.Array.Accelerate.Math.FFT.CUDA
  where

import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Twine

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Data.Complex

import Data.Array.Accelerate.CUDA.Foreign
import Data.Array.Accelerate.Array.Sugar                          as S hiding ( allocateArray )
import Data.Array.Accelerate.Type

import Control.Applicative
import Control.Concurrent.MVar
import Control.Exception
import Data.Maybe
import Foreign.CUDA.FFT
import Foreign.Storable
import System.IO.Unsafe
import System.Mem.Weak
import Text.Printf
import qualified Foreign.CUDA.Analysis                            as CUDA
import qualified Foreign.CUDA.Types                               as CUDA
import qualified Foreign.CUDA.Driver                              as CUDA
import qualified Foreign.CUDA.Driver.Context                      as Context
import Prelude                                                    as P


-- FFT using the CUFFT library to enable high performance for the CUDA backend of
-- Accelerate. The implementation works on all arrays of rank less than or equal
-- to 3. The result is un-normalised.
--
fft :: forall e sh. (Shape sh, Elt e, IsFloating e)
    => Mode
    -> sh
    -> CUDAForeignAcc (Array sh (Complex e)) (Array sh (Complex e))
fft mode sh = CUDAForeignAcc name foreignFFT
  where
    -- Plan the FFT.
    --
    -- Use 'unsafePerformIO' so that this is not performed on every invocation,
    -- if this operation is embedded in a reusable AST (e.g. run1). Note that
    -- FFT plans do not appear to be associated with a specific device context.
    --
    hndl = unsafePerformIO $ do
      plan <- case shapeToList sh of
                [width]                -> plan1D              width transform 1
                [width, height]        -> plan2D       height width transform
                [width, height, depth] -> plan3D depth height width transform
                _                      -> error "accelerate-fft cannot use CUFFT for arrays of dimensions higher than 3"
      addFinalizer plan (destroy plan)
      return plan

    -- Retrieve the CUDA functions to convert between SoA and AoS
    -- representations. These are specific to a given execution context. We use
    -- 'unsafePerformIO' again to avoid a bit of work if this operation gets
    -- embedded into the AST.
    --
    (a2c, c2a)
      | FloatingDict <- floatingDict (floatingType :: FloatingType e)
      = unsafePerformIO $ do
          let sz = sizeOf (undefined::e)
          ctx <- fromMaybe (error "could not determine current CUDA context") `fmap` Context.get
          mdl <- modifyMVar _ptx_twine_mdl $ \ms -> do
                   case lookup (ctx, sz) ms of
                     Just m  -> return (ms, m)
                     Nothing -> do
                       m <- CUDA.loadData $ case sz of
                                              4 -> ptx_twine_f32
                                              8 -> ptx_twine_f64
                                              _ -> error "I don't know what architecture I am"
                       return (((ctx,sz), m) : ms, m)
          --
          (,) <$> CUDA.getFun mdl "interleave"
              <*> CUDA.getFun mdl "deinterleave"

    -- Execute the FFT on the device, including marshalling data between AoS and
    -- SoA representations.
    --
    foreignFFT :: CUDA.Stream -> Array sh (Complex e) -> CIO (Array sh (Complex e))
    foreignFFT stream ain
      | FloatingDict  <- floatingDict (floatingType :: FloatingType e)
      = do
          let n  = S.size (S.shape ain)
          --
          ain_tmp  <- allocateArray (Z :. 2*n) :: CIO (Vector e)  -- these are really AoS (Vec2 Float) type
          aout_tmp <- allocateArray (Z :. 2*n) :: CIO (Vector e)
          aout     <- allocateArray sh
          --
          withComplexArray ain      stream $ \din_re  din_im  -> do
          withComplexArray aout     stream $ \dout_re dout_im -> do
          withScalarArray  ain_tmp  stream $ \din_tmp         -> do
          withScalarArray  aout_tmp stream $ \dout_tmp        -> do
            liftIO $ do
              dev         <- Context.device
              prp         <- CUDA.props dev
              regs        <- CUDA.requires a2c CUDA.NumRegs     -- assume same for c2a
              let
                  blockSize = 256
                  sharedMem = 0
                  maxBlocks = CUDA.maxResidentBlocks prp blockSize regs sharedMem
                  numBlocks = maxBlocks `P.min` ((n + blockSize - 1) `div` blockSize)

              -- Marshall data into AoS format
              CUDA.launchKernel a2c (numBlocks,1,1) (blockSize,1,1) sharedMem (Just stream) [CUDA.VArg din_tmp, CUDA.VArg din_re, CUDA.VArg din_im, CUDA.IArg (P.fromIntegral n)]

              -- Execute the FFT (out-of-place)
              setStream hndl stream
              executeFFT din_tmp dout_tmp

              -- Convert data into SoA format
              CUDA.launchKernel c2a (numBlocks,1,1) (blockSize,1,1) sharedMem (Just stream) [CUDA.VArg dout_re, CUDA.VArg dout_im, CUDA.VArg dout_tmp, CUDA.IArg (P.fromIntegral n)]
              return aout


    executeFFT :: CUDA.DevicePtr e -> CUDA.DevicePtr e -> IO ()
    executeFFT iptr optr
      = case (floatingType :: FloatingType e) of
          TypeFloat{}   -> execC2C hndl iptr optr sign
          TypeDouble{}  -> execZ2Z hndl iptr optr sign
          TypeCFloat{}  -> execC2C hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign
          TypeCDouble{} -> execZ2Z hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign

    withComplexArray
        :: forall sh' a.
           Array sh' (Complex e)
        -> CUDA.Stream
        -> (CUDA.DevicePtr e -> CUDA.DevicePtr e -> CIO a)
        -> CIO a
    withComplexArray (Array sh' adata) s k
      | AD_Pair (AD_Pair AD_Unit a1) a2 <- adata
      = withScalarArray (Array sh' a1 :: Array sh' e) s $ \d1 ->
        withScalarArray (Array sh' a2 :: Array sh' e) s $ \d2 ->
          k d1 d2
      --
      | otherwise
      = error "Always code as if the person who ends up maintaining your code is a violent psychopath who knows where you live."
          -- John Woods

    withScalarArray
        :: Array sh' e
        -> CUDA.Stream
        -> (CUDA.DevicePtr e -> CIO a)
        -> CIO a
    withScalarArray v s k
      = case (floatingType :: FloatingType e) of
          TypeFloat{}   -> withDevicePtr v s k
          TypeDouble{}  -> withDevicePtr v s k
          TypeCFloat{}  -> withDevicePtr v s (k . CUDA.castDevPtr)
          TypeCDouble{} -> withDevicePtr v s (k . CUDA.castDevPtr)

    withDevicePtr
        :: DevicePtrs (EltRepr e) ~ CUDA.DevicePtr b
        => Array sh' e
        -> CUDA.Stream
        -> (CUDA.DevicePtr b -> CIO a)
        -> CIO a
    withDevicePtr v s = withDevicePtrs v (Just s)

    sign :: Int
    sign = signOfMode mode

    name :: String
    name = printf "cufftExec%s.DIM%d" (show transform) (rank sh)

    transform = case (floatingType :: FloatingType e) of
                  TypeFloat{}   -> C2C
                  TypeDouble{}  -> Z2Z
                  TypeCFloat{}  -> C2C
                  TypeCDouble{} -> Z2Z


-- Functions for converting between AoS and SoA representation of complex
-- numbers, which is specific to a given execution context and type (float vs.
-- double).
--
-- Because we key this on a Context which is created as necessary via 'get', we
-- aren't keeping the context alive unnecessarily, but by the same token we have
-- no way to determine when a context becomes dead and thus remove entries from
-- this table which are no longer valid.
--
_ptx_twine_mdl :: MVar [((CUDA.Context, Int), CUDA.Module)]
_ptx_twine_mdl = unsafePerformIO $ do
  mv <- newMVar []
  _  <- mkWeakMVar mv
      $ withMVar mv
      $ mapM_ (\((ctx,_), mdl) -> bracket_ (CUDA.push ctx) CUDA.pop (CUDA.unload mdl))  -- what if the context is dead?
  return mv

