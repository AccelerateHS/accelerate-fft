{-# LANGUAGE GADTs               #-}
{-# LANGUAGE MagicHash           #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.Native.Base
-- Copyright   : [2017..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.Native.Base
  where

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Array.Unique
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.Representation.Array
import Data.Array.Accelerate.Representation.Shape
import Data.Array.Accelerate.Sugar.Elt
import Data.Primitive.Vec

import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Type
import Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix

import Data.Array.CArray.Base                                       ( CArray(..) )
import Math.FFT.Base                                                ( Sign(..), Flag, measure, preserveInput )

import Data.Bits
import Foreign.ForeignPtr
import Text.Printf
import Prelude                                                      as P


signOf :: Mode -> Sign
signOf Forward = DFTForward
signOf _       = DFTBackward

flags :: Flag
flags = measure .|. preserveInput

nameOf :: forall sh. Mode -> ShapeR sh -> String
nameOf Forward shR = printf "FFTW.dft%dD"  (rank shR)
nameOf _       shR = printf "FFTW.idft%dD" (rank shR)


-- /O(1)/ Convert a CArray to an Accelerate array
--
{-# INLINE fromCArray #-}
fromCArray
    :: forall ix sh e. (IxShapeR (EltR ix) ~ sh, Elt ix)
    => ShapeR sh
    -> NumericR e
    -> CArray ix (Complex e)
    -> IO (Array sh (Vec2 e))
fromCArray shR eR (CArray lo hi _ fp) = do
  --
  sh <- return $ rangeToShape shR (toIxShapeR lo, toIxShapeR hi) :: IO sh
  ua <- newUniqueArray (castForeignPtr fp :: ForeignPtr e)
  --
  case eR of
    NumericRfloat32 -> return $ Array sh ua
    NumericRfloat64 -> return $ Array sh ua

-- /O(1)/ Use an Accelerate array as a CArray
--
{-# INLINE withCArray #-}
withCArray
    :: forall ix sh e a. (IxShapeR (EltR ix) ~ sh, Elt ix)
    => ShapeR sh
    -> NumericR e
    -> Array sh (Vec2 e)
    -> (CArray ix (Complex e) -> IO a)
    -> IO a
withCArray shR eR (Array sh adata) f =
  let
      (lo, hi)  = shapeToRange shR sh
      wrap fp   = CArray (fromIxShapeR lo) (fromIxShapeR hi) (size shR sh) (castForeignPtr fp)
  in
  withArrayData eR adata (f . wrap)


-- Use underlying array pointers
--
{-# INLINE withArrayData #-}
withArrayData
    :: NumericR e
    -> ArrayData (Vec2 e)
    -> (ForeignPtr e -> IO a)
    -> IO a
withArrayData NumericRfloat32 ua = withLifetime (uniqueArrayData ua)
withArrayData NumericRfloat64 ua = withLifetime (uniqueArrayData ua)

