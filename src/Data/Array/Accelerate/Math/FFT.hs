{-# LANGUAGE CPP                 #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE EmptyDataDecls      #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT
-- Copyright   : [2012..2017] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
--               [2013..2017] Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- For performance, compile against the foreign library bindings (using any
-- number of '-fllvm-ptx', and '-fllvm-cpu' for the accelerate-llvm-ptx, and
-- accelerate-llvm-native backends, respectively).
--

module Data.Array.Accelerate.Math.FFT (

  Mode(..),
  Numeric,
  fft,

  fft1D,
  fft2D,
  fft3D,

) where

import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT.Type
import Data.Array.Accelerate.Math.FFT.Mode
import qualified Data.Array.Accelerate.Sugar.Shape                  as A ( rank, shapeR )
import qualified Data.Array.Accelerate.Math.FFT.Adhoc               as Adhoc

#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
import qualified Data.Array.Accelerate.Math.FFT.LLVM.Native         as Native
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.Math.FFT.LLVM.PTX            as PTX
#endif

import Prelude                                                      as P


-- | Discrete Fourier Transform along the innermost dimension of an array.
--
-- Notes for FFI implementations:
--
--   * fftw supports arrays of dimension 1-5
--   * cuFFT supports arrays of dimension 1-3
--
-- The pure implementation will be used otherwise.
--
fft :: forall sh e. (Shape sh, Slice sh, Numeric e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fft mode arr
  = let
        scale = A.fromIntegral (indexHead (shape arr))
        rank  = A.rank @(sh:.Int)
        shR   = A.shapeR @(sh:.Int)
        eR    = numericR @e
        go    =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  (if rank P.<= 5 then foreignAcc (Native.fft mode shR eR) else id) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  (if rank P.<= 3 then foreignAcc (PTX.fft    mode shR eR) else id) $
#endif
                  Adhoc.fft mode
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr


-- Vector Transform
-- ----------------

-- | Discrete Fourier Transform of a vector.
--
fft1D :: forall e. Numeric e
      => Mode
      -> Acc (Array DIM1 (Complex e))
      -> Acc (Array DIM1 (Complex e))
fft1D mode arr
  = let
        scale   = A.fromIntegral (A.length arr)
        eR      = numericR @e
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft1D mode eR) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft1D    mode eR) $
#endif
                  Adhoc.fft mode
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr


-- Matrix Transform
-- ----------------

-- | Discrete Fourier Transform of a matrix.
--
fft2D :: forall e. Numeric e
      => Mode
      -> Acc (Array DIM2 (Complex e))
      -> Acc (Array DIM2 (Complex e))
fft2D mode arr
  = let
        scale   = A.fromIntegral (A.size arr)
        eR      = numericR @e
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft2D mode eR) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft2D    mode eR) $
#endif
                  fft'

        fft' a  = A.transpose . Adhoc.fft mode
              >-> A.transpose . Adhoc.fft mode
                $ a
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr


-- Cube Transform
-- --------------

-- | Discrete Fourier Transform of a 3D array.
--
fft3D :: forall e. Numeric e
      => Mode
      -> Acc (Array DIM3 (Complex e))
      -> Acc (Array DIM3 (Complex e))
fft3D mode arr
  = let scale   = A.fromIntegral (A.size arr)
        eR      = numericR @e
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft3D mode eR) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft3D    mode eR) $
#endif
                  fft'

        fft' a  = rotate3D . Adhoc.fft mode
              >-> rotate3D . Adhoc.fft mode
              >-> rotate3D . Adhoc.fft mode
                $ a
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr


rotate3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
rotate3D arr = backpermute sh rot arr
  where
    sh :: Exp DIM3
    sh =
      let Z :. z :. y :. x = unlift (shape arr) :: Z :. Exp Int :. Exp Int :. Exp Int
      in  index3 y x z
    --
    rot :: Exp DIM3 -> Exp DIM3
    rot ix =
      let Z :. z :. y :. x = unlift ix          :: Z :. Exp Int :. Exp Int :. Exp Int
      in  index3 x z y

{--
-- Rank-generalised Cooley-Tuckey DFT
--
-- We require the innermost dimension be passed as a Haskell value because we
-- can't do divide-and-conquer recursion directly in the meta-language.
--
fft :: forall sh e. (Slice sh, Shape sh, A.RealFloat e, A.FromIntegral Int e)
    => e
    -> sh
    -> Int
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fft sign sh sz arr
  | P.any (P.not . isPow2) (shapeToList (sh:.sz))
  = error $ printf "fft: array dimensions must be powers-of-two, but are: %s" (showShape (sh:.sz))
  --
  | otherwise
  = go sz 0 1
  where
    go :: Int -> Int -> Int -> Acc (Array (sh:.Int) (Complex e))
    go len offset stride
      | len P.== 2
      = A.generate (constant (sh :. len)) swivel

      | otherwise
      = combine
          (go (len `div` 2) offset            (stride * 2))
          (go (len `div` 2) (offset + stride) (stride * 2))

      where
        len'    = the (unit (constant len))
        offset' = the (unit (constant offset))
        stride' = the (unit (constant stride))

        swivel ix =
          let sh' :. sz' = unlift ix :: Exp sh :. Exp Int
          in
          sz' A.== 0 ? ( (arr ! lift (sh' :. offset')) + (arr ! lift (sh' :. offset' + stride'))
          {-  A.== 1-} , (arr ! lift (sh' :. offset')) - (arr ! lift (sh' :. offset' + stride')) )

        combine evens odds =
          let odds' = A.generate (A.shape odds) (\ix -> twiddle len' (indexHead ix) * odds!ix)
          in
          append (A.zipWith (+) evens odds') (A.zipWith (-) evens odds')

        twiddle n' i' =
          let n = A.fromIntegral n'
              i = A.fromIntegral i'
              k = 2*pi*i/n
          in
          lift ( cos k :+ A.constant sign * sin k )


-- Append two arrays. This is a specialised version of (A.++) which does not do
-- bounds checking or intersection.
--
append
    :: forall sh e. (Slice sh, Shape sh, Elt e)
    => Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
append xs ys
  = let sh :. n = unlift (A.shape xs)     :: Exp sh :. Exp Int
        _  :. m = unlift (A.shape ys)     :: Exp sh :. Exp Int
    in
    generate (lift (sh :. n+m))
             (\ix -> let sz :. i = unlift ix :: Exp sh :. Exp Int
                     in  i A.< n ? (xs ! lift (sz:.i), ys ! lift (sz:.i-n) ))

isPow2 :: Int -> Bool
isPow2 0 = True
isPow2 1 = False
isPow2 x = x .&. (x-1) P.== 0
--}

