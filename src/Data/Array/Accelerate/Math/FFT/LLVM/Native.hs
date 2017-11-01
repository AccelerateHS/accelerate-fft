{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE PatternGuards       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.Native
-- Copyright   : [2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.Native (

  fft,
  fft1D,
  fft2D,
  fft3D,

) where

import Data.Array.Accelerate.Math.FFT.Mode

import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Analysis.Match                         as A
import Data.Array.Accelerate.Array.Data                             as A
import Data.Array.Accelerate.Array.Sugar                            as S
import Data.Array.Accelerate.Array.Unique                           as A
import Data.Array.Accelerate.Data.Complex                           as A
import Data.Array.Accelerate.Error                                  as A
import Data.Array.Accelerate.Type                                   as A

import Data.Array.Accelerate.LLVM.Native.Foreign

import Data.Ix                                                      ( Ix )
import Data.Array.CArray                                            ( CArray )
import qualified Data.Array.CArray                                  as C

import Math.FFT.Base                                                ( FFTWReal, Sign(..), Flag, measure, destroyInput )
import qualified Math.FFT                                           as FFT

import Foreign.Ptr
import Foreign.Storable
import Foreign.Storable.Complex                                     ()

import Data.Bits
import Data.Typeable
import Text.Printf
import Prelude                                                      as P


fft :: forall sh e. (Shape sh, Elt e, IsFloating e)
    => Mode
    -> ForeignAcc (Array sh (Complex e) -> Array sh (Complex e))
fft mode
  = ForeignAcc (nameOf mode (undefined::sh))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . go
      TypeDouble{}  -> liftIO . go
      TypeCFloat{}  -> liftIO . go
      TypeCDouble{} -> liftIO . go
  where
    go :: (ArrayPtrs e ~ Ptr r, FFTWReal r) => Array sh (Complex e) -> IO (Array sh (Complex e))
    go | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM1) = liftAtoC (FFT.dftGU (signOf mode) flags [0] `ix` (undefined :: Int))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM2) = liftAtoC (FFT.dftGU (signOf mode) flags [1] `ix` (undefined :: (Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM3) = liftAtoC (FFT.dftGU (signOf mode) flags [2] `ix` (undefined :: (Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM4) = liftAtoC (FFT.dftGU (signOf mode) flags [3] `ix` (undefined :: (Int,Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM5) = liftAtoC (FFT.dftGU (signOf mode) flags [4] `ix` (undefined :: (Int,Int,Int,Int,Int)))
       | otherwise = $internalError "fft" "only for 1D..5D inner-dimension transforms"
    --
    ix :: (a i r -> a i r) -> i -> (a i r -> a i r)
    ix f _ = f


fft1D :: forall e. (Elt e, IsFloating e)
      => Mode
      -> ForeignAcc (Vector (Complex e) -> Vector (Complex e))
fft1D mode
  = ForeignAcc (nameOf mode (undefined::DIM1))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . liftAtoC go
      TypeDouble{}  -> liftIO . liftAtoC go
      TypeCFloat{}  -> liftIO . liftAtoC go
      TypeCDouble{} -> liftIO . liftAtoC go
  where
    go :: FFTWReal r => CArray Int (Complex r) -> CArray Int (Complex r)
    go = FFT.dftGU (signOf mode) flags [0]


fft2D :: forall e. (Elt e, IsFloating e)
      => Mode
      -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
fft2D mode
  = ForeignAcc (nameOf mode (undefined::DIM2))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . liftAtoC go
      TypeDouble{}  -> liftIO . liftAtoC go
      TypeCFloat{}  -> liftIO . liftAtoC go
      TypeCDouble{} -> liftIO . liftAtoC go
  where
    go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1]


fft3D :: forall e. (Elt e, IsFloating e)
      => Mode
      -> ForeignAcc (Array DIM3 (Complex e) -> Array DIM3 (Complex e))
fft3D mode
  = ForeignAcc (nameOf mode (undefined::DIM3))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . liftAtoC go
      TypeDouble{}  -> liftIO . liftAtoC go
      TypeCFloat{}  -> liftIO . liftAtoC go
      TypeCDouble{} -> liftIO . liftAtoC go
  where
    go :: FFTWReal r => CArray (Int,Int,Int) (Complex r) -> CArray (Int,Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1,2]


signOf :: Mode -> Sign
signOf Forward = DFTForward
signOf _       = DFTBackward

flags :: Flag
flags = measure .|. destroyInput

nameOf :: forall sh. Shape sh => Mode -> sh -> String
nameOf Forward _ = printf "FFTW.dft%dD"  (rank (undefined::sh))
nameOf _       _ = printf "FFTW.idft%dD" (rank (undefined::sh))


-- | Lift an operation on CArray into an operation on Accelerate arrays
--
liftAtoC
    :: (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Ix ix, Elt ix, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => (CArray ix (Complex e') -> CArray ix (Complex e'))
    -> Array sh (Complex e)
    -> IO (Array sh (Complex e))
liftAtoC f a = c2a . f =<< a2c a


-- | Convert a multidimensional Accelerate array of complex numbers into
-- a packed CArray of complex numbers suitable for use by FFTW.
--
a2c :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => Array sh (Complex e)
    -> IO (CArray ix (Complex e'))
a2c arr
  | FloatingDict <- floatingDict (floatingType :: FloatingType e)
  = let
        (lo,hi) = shapeToRange (arrayShape arr)
        bnds    = (fromIxShapeRepr lo, fromIxShapeRepr hi)
        n       = S.size (arrayShape arr)
    in
    C.createCArray       bnds $ \p_cs      ->
    withComplexArrayPtrs arr  $ \p_re p_im ->
      let
          -- TLM: Should we execute this in parallel using the worker threads of
          -- the current target? (Native)
          go !i | i P.>= n = return ()
          go !i            = do
            re <- peekElemOff p_re i
            im <- peekElemOff p_im i
            pokeElemOff p_cs i (re :+ im)
            go (i+1)
      in
      go 0


-- | Convert a packed CArray of complex numbers into an unzipped (SoA)
-- multidimensional Accelerate array of complex numbers.
--
c2a :: forall ix sh e e'. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Ix ix, Elt ix, Shape sh, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => CArray ix (Complex e')
    -> IO (Array sh (Complex e))
c2a carr
  | FloatingDict <- floatingDict (floatingType :: FloatingType e)
  = let
        (lo,hi) = C.bounds carr
        n       = C.rangeSize (lo,hi)
        sh      = rangeToShape (toIxShapeRepr lo, toIxShapeRepr hi)
    in do
      arr <- allocateArray sh
      C.withCArray carr        $ \p_cs      -> do
      withComplexArrayPtrs arr $ \p_re p_im -> do
        let
            -- TLM: Should we execute this in parallel using the worker threads
            -- of the current target? (Native)
            go !i | i P.>= n = return ()
            go !i            = do
              re :+ im <- peekElemOff p_cs i
              pokeElemOff p_re i re
              pokeElemOff p_im i im
              go (i+1)
        --
        go 0
        return arr


-- Converting between Accelerate multidimensional shapes/indices and those used
-- by the CArray package (Data.Ix)
--

type family IxShapeRepr e where
  IxShapeRepr ()    = ()
  IxShapeRepr Int   = ((),Int)
  IxShapeRepr (t,h) = (IxShapeRepr t, h)

fromIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => sh
    -> ix
fromIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> IxShapeRepr ix' -> ix'
    go UnitTuple                                                 ()     = ()
    go (PairTuple tt _)                                          (t, h) = (go tt t, h)
    go (SingleTuple (NumScalarType (IntegralNumType TypeInt{}))) ((),h) = h
    go _ _
      = $internalError "fromIxShapeRepr" "expected Int dimensions"

toIxShapeRepr
    :: forall ix sh. (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Elt ix)
    => ix
    -> sh
toIxShapeRepr = liftToElt (go (eltType (undefined::ix)))
  where
    go :: forall ix'. TupleType ix' -> ix' -> IxShapeRepr ix'
    go UnitTuple        ()                                             = ()
    go (SingleTuple     (NumScalarType (IntegralNumType TypeInt{}))) h = ((), h)
    go (PairTuple tt _) (t, h)                                         = (go tt t, h)
    go _ _
      = error "toIxShapeRepr: not a valid Data.Ix index"


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

-- withScalarArrayPtrs
--     :: forall sh e a. IsFloating e
--     => Array sh e
--     -> (ArrayPtrs e -> IO a)
--     -> IO a
-- withScalarArrayPtrs (Array _ adata) =
--   case floatingType :: FloatingType e of
--     TypeFloat{}   -> withArrayData arrayElt adata
--     TypeDouble{}  -> withArrayData arrayElt adata
--     TypeCFloat{}  -> withArrayData arrayElt adata
--     TypeCDouble{} -> withArrayData arrayElt adata

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

