{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
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

import Data.Array.Accelerate.Analysis.Match
import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.Type

import Data.Array.Accelerate.LLVM.PTX.Foreign

import Foreign.CUDA.Ptr                                             ( DevicePtr )

import Data.Typeable
import Foreign.Ptr
import Foreign.Storable


type family DevicePtrs e :: *

type instance DevicePtrs ()      = ()
type instance DevicePtrs Float   = DevicePtr Float
type instance DevicePtrs Double  = DevicePtr Double
type instance DevicePtrs CFloat  = DevicePtr Float
type instance DevicePtrs CDouble = DevicePtr Double
type instance DevicePtrs (a,b)   = (DevicePtrs a, DevicePtrs b)


-- | Dig out the two device pointers for an unzipped array of complex numbers.
--
withComplexArrayPtrs
    :: forall sh e a. IsFloating e
    => Array sh (Complex e)
    -> Stream
    -> (DevicePtrs e -> DevicePtrs e -> LLVM PTX a)
    -> LLVM PTX a
withComplexArrayPtrs (Array _ adata) st k
  | AD_Pair (AD_Pair AD_Unit ad1) ad2 <- adata
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withArrayData arrayElt ad1 st $ \p1 -> withArrayData arrayElt ad2 st $ \p2 -> k p1 p2
      TypeDouble{}  -> withArrayData arrayElt ad1 st $ \p1 -> withArrayData arrayElt ad2 st $ \p2 -> k p1 p2
      TypeCDouble{} -> withArrayData arrayElt ad1 st $ \p1 -> withArrayData arrayElt ad2 st $ \p2 -> k p1 p2
      TypeCFloat{}  -> withArrayData arrayElt ad1 st $ \p1 -> withArrayData arrayElt ad2 st $ \p2 -> k p1 p2

-- | Dig out the device pointer for a scalar array
--
withScalarArrayPtr
    :: forall sh e a. IsFloating e
    => Array sh e
    -> Stream
    -> (DevicePtrs e -> LLVM PTX a)
    -> LLVM PTX a
withScalarArrayPtr (Array _ ad) st k
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withArrayData arrayElt ad st $ \p -> k p
      TypeDouble{}  -> withArrayData arrayElt ad st $ \p -> k p
      TypeCDouble{} -> withArrayData arrayElt ad st $ \p -> k p
      TypeCFloat{}  -> withArrayData arrayElt ad st $ \p -> k p

withArrayData
    :: (Typeable e, Typeable a, ArrayElt e, Storable a, ArrayPtrs e ~ Ptr a)
    => ArrayEltR e
    -> ArrayData e
    -> Stream
    -> (DevicePtr a -> LLVM PTX b)
    -> LLVM PTX b
withArrayData _ ad s k =
  withDevicePtr ad $ \p -> do
    r <- k p
    e <- checkpoint s
    return (Just e,r)


withLifetime' :: Lifetime a -> (a -> LLVM PTX b) -> LLVM PTX b
withLifetime' l k = do
  r <- k (unsafeGetValue l)
  liftIO $ touchLifetime l
  return r


-- Match reified shape types
--
matchShapeType
    :: forall sh sh'. (Shape sh, Shape sh')
    => sh
    -> sh'
    -> Maybe (sh :~: sh')
matchShapeType _ _
  | Just Refl <- matchTupleType (eltType (undefined::sh)) (eltType (undefined::sh'))
  = gcast Refl

matchShapeType _ _
  = Nothing

