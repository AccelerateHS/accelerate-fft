{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
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

import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Type


-- Converting between Accelerate multidimensional shapes/indices and those used
-- by the CArray package (Data.Ix)
--
type family IxShapeRepr e where
  IxShapeRepr ()    = ()
  IxShapeRepr Int   = ((),Int)
  IxShapeRepr (t,h) = (IxShapeRepr t, h)

{-# INLINE fromIxShapeRepr #-}
fromIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => sh
    -> ix
fromIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> IxShapeRepr ix' -> ix'
    go TypeRunit                                                                    ()     = ()
    go (TypeRpair tt _)                                                             (t, h) = (go tt t, h)
    go (TypeRscalar (SingleScalarType (NumSingleType (IntegralNumType TypeInt{})))) ((),h) = h
    go _ _
      = $internalError "fromIxShapeRepr" "expected Int dimensions"

{-# INLINE toIxShapeRepr #-}
toIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => ix
    -> sh
toIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> ix' -> IxShapeRepr ix'
    go TypeRunit        ()                                                                = ()
    go (TypeRscalar     (SingleScalarType (NumSingleType (IntegralNumType TypeInt{})))) h = ((), h)
    go (TypeRpair tt _) (t, h)                                                            = (go tt t, h)
    go _ _
      = $internalError "toIxShapeRepr" "not a valid Data.Ix index"

