{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeFamilies        #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix
  where

import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Representation.Type
import Data.Array.Accelerate.Sugar.Elt
import Data.Array.Accelerate.Type


-- Converting between Accelerate multidimensional shapes/indices and those used
-- by the CArray package (Data.Ix)
--
type family IxShapeR e where
  IxShapeR ()    = ()
  IxShapeR Int   = ((),Int)
  IxShapeR (t,h) = (IxShapeR t, h)

{-# INLINE fromIxShapeR #-}
fromIxShapeR
    :: forall ix sh. (HasCallStack, IxShapeR (EltR ix) ~ sh, Elt ix)
    => sh
    -> ix
fromIxShapeR = toElt . go (eltR @ix)
  where
    go :: forall ix'. TypeR ix' -> IxShapeR ix' -> ix'
    go TupRunit                                                                    ()     = ()
    go (TupRpair tt _)                                                             (t, h) = (go tt t, h)
    go (TupRsingle (SingleScalarType (NumSingleType (IntegralNumType TypeInt{})))) ((),h) = h
    go _ _
      = internalError "expected Int dimensions"

{-# INLINE toIxShapeR #-}
toIxShapeR
    :: forall ix sh. (HasCallStack, IxShapeR (EltR ix) ~ sh, Elt ix)
    => ix
    -> sh
toIxShapeR = go (eltR @ix) . fromElt
  where
    go :: forall ix'. TypeR ix' -> ix' -> IxShapeR ix'
    go TupRunit        ()                                                                = ()
    go (TupRsingle     (SingleScalarType (NumSingleType (IntegralNumType TypeInt{})))) h = ((), h)
    go (TupRpair tt _) (t, h)                                                            = (go tt t, h)
    go _ _
      = internalError "not a valid Data.Ix index"

