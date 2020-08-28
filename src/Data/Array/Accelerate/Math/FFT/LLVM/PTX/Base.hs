{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base
-- Copyright   : [2017..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base
  where

import Data.Array.Accelerate.Math.FFT.Type

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.Representation.Array
import Data.Array.Accelerate.Type
import Data.Primitive.Vec

import Data.Array.Accelerate.LLVM.PTX.Foreign

import Foreign.CUDA.Ptr                                             ( DevicePtr )


{-# INLINE withArray #-}
withArray
    :: NumericR e
    -> Array sh (Vec2 e)
    -> Stream
    -> (DevicePtr e -> LLVM PTX b)
    -> LLVM PTX b
withArray eR (Array _ adata) = withArrayData eR adata

{-# INLINE withArrayData #-}
withArrayData
    :: NumericR e
    -> ArrayData (Vec2 e)
    -> Stream
    -> (DevicePtr e -> LLVM PTX b)
    -> LLVM PTX b
withArrayData NumericRfloat32 ad s k =
  withDevicePtr (singleType @Float) ad $ \p -> do
    r <- k p
    e <- waypoint s
    return (Just e,r)
withArrayData NumericRfloat64 ad s k =
  withDevicePtr (singleType @Double) ad $ \p -> do
    r <- k p
    e <- waypoint s
    return (Just e, r)

{-# INLINE withLifetime' #-}
withLifetime' :: Lifetime a -> (a -> LLVM PTX b) -> LLVM PTX b
withLifetime' l k = do
  r <- k (unsafeGetValue l)
  liftIO $ touchLifetime l
  return r

