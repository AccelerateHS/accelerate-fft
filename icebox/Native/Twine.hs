{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE GADTs                    #-}
{-# LANGUAGE PatternGuards            #-}
{-# LANGUAGE ScopedTypeVariables      #-}
{-# LANGUAGE TemplateHaskell          #-}
{-# LANGUAGE TypeFamilies             #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.Native.Twine
-- Copyright   : [2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.Native.Twine (

  interleave,
  deinterleave,

) where

import Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix

import Data.Array.Accelerate.Array.Data
import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Array.Unique
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Error
import Data.Array.Accelerate.Type

import Foreign.Ptr
import Foreign.Storable

import Data.Ix                                                      ( Ix )
import Data.Array.CArray                                            ( CArray )
import qualified Data.Array.CArray                                  as C


-- | Convert a multidimensional Accelerate array of complex numbers into
-- a packed CArray of complex numbers suitable for use by FFTW.
--
-- TODO: should we execute this in parallel?
--
interleave
    :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => Array sh (Complex e)
    -> IO (CArray ix (Complex e'))
interleave arr =
  let
      (lo,hi) = shapeToRange (shape arr)
      bnds    = (fromIxShapeRepr lo, fromIxShapeRepr hi)
      n       = size (shape arr)
  in
  C.createCArray       bnds $ \p_cs      ->
  withComplexArrayPtrs arr  $ \p_re p_im ->
    case floatingType :: FloatingType e of
      TypeFloat{}   -> c_interleave_f32 0 n p_cs p_re p_im
      TypeDouble{}  -> c_interleave_f64 0 n p_cs p_re p_im
      TypeCFloat{}  -> c_interleave_f32 0 n p_cs p_re p_im
      TypeCDouble{} -> c_interleave_f64 0 n p_cs p_re p_im


-- | Convert a packed CArray of complex numbers into an unzipped (SoA)
-- multidimensional Accelerate array of complex numbers.
--
-- TODO: should we execute this in parallel?
--
deinterleave
    :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => CArray ix (Complex e')
    -> IO (Array sh (Complex e))
deinterleave carr = do
  let
      (lo,hi) = C.bounds carr
      n       = C.rangeSize (lo,hi)
      sh      = rangeToShape (toIxShapeRepr lo, toIxShapeRepr hi)
  --
  arr <- allocateArray sh
  C.withCArray carr        $ \p_cs      -> do
  withComplexArrayPtrs arr $ \p_re p_im -> do
    () <- case floatingType :: FloatingType e of
            TypeFloat{}   -> c_deinterleave_f32 0 n p_re p_im p_cs
            TypeDouble{}  -> c_deinterleave_f64 0 n p_re p_im p_cs
            TypeCFloat{}  -> c_deinterleave_f32 0 n p_re p_im p_cs
            TypeCDouble{} -> c_deinterleave_f64 0 n p_re p_im p_cs
    --
    return arr


-- Dig out the underlying pointers of the Accelerate SoA data structure
--
withComplexArrayPtrs
    :: forall sh e a. IsFloating e
    => Array sh (Complex e)
    -> (ArrayPtrs e -> ArrayPtrs e -> IO a)
    -> IO a
withComplexArrayPtrs (Array _ adata) k
  | AD_Pair (AD_Pair AD_Unit ad1) ad2 <- adata
  = case floatingType :: FloatingType e of
      TypeFloat{}   -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeDouble{}  -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeCFloat{}  -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2
      TypeCDouble{} -> withArrayData arrayElt ad1 $ \p1 -> withArrayData arrayElt ad2 $ \p2 -> k p1 p2

withArrayData
    :: (ArrayPtrs e ~ Ptr a)
    => ArrayEltR e
    -> ArrayData e
    -> (Ptr a -> IO b)
    -> IO b
withArrayData ArrayEltRfloat   (AD_Float   ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRdouble  (AD_Double  ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRcfloat  (AD_CFloat  ua) = withUniqueArrayPtr ua
withArrayData ArrayEltRcdouble (AD_CDouble ua) = withUniqueArrayPtr ua
withArrayData _ _ =
  $internalError "withArrayData" "expected array of [C]Float or [C]Double"


foreign import ccall unsafe "interleave_f32"
  c_interleave_f32 :: Int -> Int -> Ptr (Complex Float) -> Ptr Float -> Ptr Float -> IO ()

foreign import ccall unsafe "interleave_f64"
  c_interleave_f64 :: Int -> Int -> Ptr (Complex Double) -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "deinterleave_f32"
  c_deinterleave_f32 :: Int -> Int -> Ptr Float -> Ptr Float -> Ptr (Complex Float) -> IO ()

foreign import ccall unsafe "deinterleave_f64"
  c_deinterleave_f64 :: Int -> Int -> Ptr Double -> Ptr Double -> Ptr (Complex Double) -> IO ()

