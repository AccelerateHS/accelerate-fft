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

import Data.Array.Accelerate                                      as A
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Twine

import Data.Array.Accelerate.CUDA.Foreign
import Data.Array.Accelerate.Array.Sugar                          as S ( shapeToList, shape, EltRepr )
import Data.Array.Accelerate.Type

import System.Mem.Weak
import System.IO.Unsafe
import Foreign.CUDA.FFT
import qualified Foreign.CUDA.Types                               as CUDA
import qualified Foreign.CUDA.Driver                              as CUDA


-- FFT using the CUFFT library to enable high performance for the CUDA backend of
-- Accelerate. The implementation works on all arrays of rank less than or equal
-- to 3. The result is un-normalised.
--
fft :: forall e sh. (Shape sh, Elt e, IsFloating e)
    => Mode
    -> sh
    -> (Acc (Array sh (Complex e)) -> Acc (Array sh (Complex e)))
    -> Acc (Array sh (Complex e))
    -> Acc (Array sh (Complex e))
fft mode sh = cudaFFT'
  where
    -- Plan the FFT.
    -- Doing this in unsafePerformIO so it is not reperformed every time the
    -- AST is evaluated.
    --
    hndl = unsafePerformIO $ do
            plan <- case shapeToList sh of
                     [width]                -> plan1D              width types 1
                     [width, height]        -> plan2D       height width types
                     [width, height, depth] -> plan3D depth height width types
                     _                      -> error "Accelerate-fft cannot use CUFFT for arrays of dimensions higher than 3"
            addFinalizer plan (destroy plan)
            return plan

    types = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> C2C
              TypeDouble{}  -> Z2Z
              TypeCFloat{}  -> C2C
              TypeCDouble{} -> Z2Z

    cudaFFT' p
      = reshape (constant sh)
      . deinterleave
      . foreignAcc ff pure
      . interleave
      . flatten
      where
        ff          = CUDAForeignAcc "foreignFFT" foreignFFT
        -- Unfortunately the pure version of the function needs to be wrapped in
        -- interleave and deinterleave to match how the foreign version works.
        --
        -- RCE: Do the interleaving and deinterleaving in foreignFFT
        --
        -- TLM: The interleaving might get fused into other parts of the
        --      computation and thus be okay. We should really support multi types
        --      such as float2 instead.
        --
        pure        = interleave . flatten . p . reshape (constant sh) . deinterleave
        sign        = signOfMode mode :: Int

        foreignFFT :: CUDA.Stream -> Array DIM1 e -> CIO (Array DIM1 e)
        foreignFFT stream arr' = do
          output <- allocateArray (S.shape arr')
          withFloatingDevicePtr arr'   stream $ \iptr -> do
          withFloatingDevicePtr output stream $ \optr -> do
            -- Execute the foreign function
            liftIO $ do
              setStream hndl stream
              execute iptr optr
              return output

        execute :: CUDA.DevicePtr e -> CUDA.DevicePtr e -> IO ()
        execute iptr optr
          = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> execC2C hndl iptr optr sign
              TypeDouble{}  -> execZ2Z hndl iptr optr sign
              TypeCFloat{}  -> execC2C hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign
              TypeCDouble{} -> execZ2Z hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign

        withFloatingDevicePtr :: Vector e -> CUDA.Stream -> (CUDA.DevicePtr e -> CIO a) -> CIO a
        withFloatingDevicePtr v s k
          = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> withSingleDevicePtr v s k
              TypeDouble{}  -> withSingleDevicePtr v s k
              TypeCFloat{}  -> withSingleDevicePtr v s (k . CUDA.castDevPtr)
              TypeCDouble{} -> withSingleDevicePtr v s (k . CUDA.castDevPtr)

        withSingleDevicePtr
            :: DevicePtrs (EltRepr e) ~ CUDA.DevicePtr b
            => Vector e
            -> CUDA.Stream
            -> (CUDA.DevicePtr b -> CIO a)
            -> CIO a
        withSingleDevicePtr v s = withDevicePtrs v (Just s)

