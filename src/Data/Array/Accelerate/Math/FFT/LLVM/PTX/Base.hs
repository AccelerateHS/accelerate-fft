{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base
  where

import Data.Array.Accelerate.Math.FFT.Type

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Lifetime

import Data.Array.Accelerate.LLVM.PTX.Foreign

import Foreign.CUDA.Ptr                                             ( DevicePtr )


{-# INLINE withArray #-}
withArray
    :: forall sh e b. Numeric e
    => Array sh (Complex e)
    -> Stream
    -> (DevicePtr e -> LLVM PTX b)
    -> LLVM PTX b
withArray (Array _ adata) = withArrayData (numericR::NumericR e) adata

{-# INLINE withArrayData #-}
withArrayData
    :: NumericR e
    -> ArrayData (EltRepr (Complex e))
    -> Stream
    -> (DevicePtr e -> LLVM PTX b)
    -> LLVM PTX b
withArrayData NumericRfloat32 (AD_Vec _ ad) s k =
  withDevicePtr ad $ \p -> do
    r <- k p
    e <- waypoint s
    return (Just e,r)
withArrayData NumericRfloat64 (AD_Vec _ ad) s k =
  withDevicePtr ad $ \p -> do
    r <- k p
    e <- waypoint s
    return (Just e, r)

{-# INLINE withLifetime' #-}
withLifetime' :: Lifetime a -> (a -> LLVM PTX b) -> LLVM PTX b
withLifetime' l k = do
  r <- k (unsafeGetValue l)
  liftIO $ touchLifetime l
  return r

