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
import Data.Array.Accelerate.Math.FFT.LLVM.Native.Ix
import Data.Array.Accelerate.Math.FFT.LLVM.Native.Twine

import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Analysis.Match                         as A
import Data.Array.Accelerate.Array.Data                             as A
import Data.Array.Accelerate.Array.Sugar                            as S
import Data.Array.Accelerate.Data.Complex                           as A
import Data.Array.Accelerate.Error                                  as A
import Data.Array.Accelerate.Type                                   as A

import Data.Array.Accelerate.LLVM.Native.Foreign

import Data.Ix                                                      ( Ix )
import Data.Array.CArray                                            ( CArray )

import Math.FFT.Base                                                ( FFTWReal, Sign(..), Flag, measure, destroyInput )
import qualified Math.FFT                                           as FFT

import Foreign.Ptr
import Foreign.Storable

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
    go | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM1) = liftCtoA (FFT.dftGU (signOf mode) flags [0] `ix` (undefined :: Int))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM2) = liftCtoA (FFT.dftGU (signOf mode) flags [1] `ix` (undefined :: (Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM3) = liftCtoA (FFT.dftGU (signOf mode) flags [2] `ix` (undefined :: (Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM4) = liftCtoA (FFT.dftGU (signOf mode) flags [3] `ix` (undefined :: (Int,Int,Int,Int)))
       | Just Refl <- matchShapeType (undefined::sh) (undefined::DIM5) = liftCtoA (FFT.dftGU (signOf mode) flags [4] `ix` (undefined :: (Int,Int,Int,Int,Int)))
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
      TypeFloat{}   -> liftIO . liftCtoA go
      TypeDouble{}  -> liftIO . liftCtoA go
      TypeCFloat{}  -> liftIO . liftCtoA go
      TypeCDouble{} -> liftIO . liftCtoA go
  where
    go :: FFTWReal r => CArray Int (Complex r) -> CArray Int (Complex r)
    go = FFT.dftGU (signOf mode) flags [0]


fft2D :: forall e. (Elt e, IsFloating e)
      => Mode
      -> ForeignAcc (Array DIM2 (Complex e) -> Array DIM2 (Complex e))
fft2D mode
  = ForeignAcc (nameOf mode (undefined::DIM2))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . liftCtoA go
      TypeDouble{}  -> liftIO . liftCtoA go
      TypeCFloat{}  -> liftIO . liftCtoA go
      TypeCDouble{} -> liftIO . liftCtoA go
  where
    go :: FFTWReal r => CArray (Int,Int) (Complex r) -> CArray (Int,Int) (Complex r)
    go = FFT.dftGU (signOf mode) flags [0,1]


fft3D :: forall e. (Elt e, IsFloating e)
      => Mode
      -> ForeignAcc (Array DIM3 (Complex e) -> Array DIM3 (Complex e))
fft3D mode
  = ForeignAcc (nameOf mode (undefined::DIM3))
  $ case floatingType :: FloatingType e of
      TypeFloat{}   -> liftIO . liftCtoA go
      TypeDouble{}  -> liftIO . liftCtoA go
      TypeCFloat{}  -> liftIO . liftCtoA go
      TypeCDouble{} -> liftIO . liftCtoA go
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
liftCtoA
    :: (IxShapeRepr (EltRepr ix) ~ EltRepr sh, Shape sh, Ix ix, Elt ix, Elt e, IsFloating e, Storable e', ArrayPtrs e ~ Ptr e')
    => (CArray ix (Complex e') -> CArray ix (Complex e'))
    -> Array sh (Complex e)
    -> IO (Array sh (Complex e))
liftCtoA f a = deinterleave . f =<< interleave a


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

