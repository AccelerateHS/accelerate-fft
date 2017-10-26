{-# LANGUAGE CPP                      #-}
{-# LANGUAGE ConstraintKinds          #-}
{-# LANGUAGE EmptyDataDecls           #-}
{-# LANGUAGE FlexibleContexts         #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE GADTs                    #-}
{-# LANGUAGE ScopedTypeVariables      #-}
{-# LANGUAGE TypeFamilies             #-}
{-# LANGUAGE TypeOperators            #-}
{-# LANGUAGE ViewPatterns             #-}
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
-- Computation of a Discrete Fourier Transform using the Cooley-Tuckey
-- algorithm. The time complexity is O(n log n) in the size of the input.
--
-- The base (default) implementation uses a naïve divide-and-conquer algorithm
-- whose absolute performance is appalling. It also requires that you know on
-- the Haskell side the size of the data being transformed, and that this is
-- a power-of-two in each dimension.
--
-- For performance, compile against the foreign library bindings (using any
-- number of '-fllvm-ptx', and '-fllvm-cpu' for the accelerate-llvm-ptx, and
-- accelerate-llvm-native backends, respectively), which have none of the above
-- restrictions.
--

module Data.Array.Accelerate.Math.FFT (

  Mode(..),
  FFTElt,
  fft1D, fft1D', fft1D_2r', fft1D_3r',
  fft2D, fft2D',
  fft3D, fft3D',
  fft

) where

import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Array.Sugar                            ( showShape, shapeToList )
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT.Mode

#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
import qualified Data.Array.Accelerate.Math.FFT.LLVM.Native         as Native
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
import qualified Data.Array.Accelerate.Math.FFT.LLVM.PTX            as PTX
#endif

import Data.Bits
import Text.Printf
import Prelude                                                      as P


-- The type of supported FFT elements; namely 'Float' and 'Double'.
--
type FFTElt e = (P.Num e, A.RealFloat e, A.FromIntegral Int e, A.IsFloating e)


-- Vector Transform
-- ----------------

-- | Discrete Fourier Transform of a vector.
--
-- The default implementation requires the array dimension to be a power of two
-- (else error).
--
fft1D :: FFTElt e
      => Mode
      -> Array DIM1 (Complex e)
      -> Acc (Array DIM1 (Complex e))
fft1D mode vec
  = fft1D' mode (arrayShape vec) (use vec)


-- | Discrete Fourier Transform of a vector.
--
-- The default implementation requires the array dimension to be a power of two.
-- The FFI-backed implementations ignore the Haskell-side size parameter (second
-- argument).
--
fft1D' :: forall e. FFTElt e
       => Mode
       -> DIM1
       -> Acc (Array DIM1 (Complex e))
       -> Acc (Array DIM1 (Complex e))
fft1D' mode (Z :. len) arr
  = let sign    = signOfMode mode :: e
        scale   = A.fromIntegral (A.length arr)
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft1D mode) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft1D mode) $
#endif
                  fft sign Z len
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr


-- Matrix Transform
-- ----------------

-- | Discrete Fourier Transform of a matrix.
--
-- The default implementation requires the array dimensions to be powers of two
-- (else error).
--
fft2D :: FFTElt e
      => Mode
      -> Array DIM2 (Complex e)
      -> Acc (Array DIM2 (Complex e))
fft2D mode arr
  = fft2D' mode (arrayShape arr) (use arr)


-- | Discrete Fourier Transform of a matrix.
--
-- The default implementation requires the array dimensions to be powers of two.
-- The FFI-backed implementations ignore the Haskell-side size parameter (second
-- argument).
--
fft2D' :: forall e. FFTElt e
       => Mode
       -> DIM2
       -> Acc (Array DIM2 (Complex e))
       -> Acc (Array DIM2 (Complex e))
fft2D' mode (Z :. height :. width) arr
  = let sign    = signOfMode mode :: e
        scale   = A.fromIntegral (A.size arr)
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft2D mode) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft2D mode) $
#endif
                  fft'

        fft' a  = A.transpose . fft sign (Z:.height) width
              >-> A.transpose . fft sign (Z:.width)  height
                $ a
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr

-- | Discrete Fourier Transform of all rows in a matrix.
--
-- The default implementation requires the row`s length to be a power of two.
-- The FFI-backed implementations ignore the Haskell-side size parameter (second
-- argument).

fft1D_2r' :: forall e sh. (FFTElt e, Shape sh, sh ~ DIM1)
         => Mode
         -> (sh :. Int)
         -> Acc (Array DIM2 (Complex e))
         -> Acc (Array DIM2 (Complex e))
fft1D_2r' mode (A.Z :. height :. width) arr
  = let sign    = signOfMode mode :: e
        (A.Z A.:. (eHeight :: A.Exp Int) A.:. (_ :: A.Exp Int))  = A.unlift $ A.shape arr
        scale   = (A.fromIntegral (A.size arr))/(A.fromIntegral $ eHeight)
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft1D_r mode) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft1D_r mode) $
#endif
                  fft sign (Z:.width) height
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr

-- | Discrete Fourier Transform of all rows in a 3D array.
--
-- The default implementation requires the row`s length to be a power of two.
-- The FFI-backed implementations ignore the Haskell-side size parameter (second
-- argument).

fft1D_3r' :: forall e sh. (FFTElt e, Shape sh, sh ~ DIM2)
         => Mode
         -> (sh :. Int)
         -> Acc (Array DIM3 (Complex e))
         -> Acc (Array DIM3 (Complex e))
fft1D_3r' mode (A.Z :. depth :. height :. width) arr
  = let sign    = signOfMode mode :: e
        (A.Z A.:. (eDepth :: A.Exp Int) A.:. (eHeight :: A.Exp Int) A.:. (_ :: A.Exp Int))  = A.unlift $ A.shape arr
        scale   = (A.fromIntegral (A.size arr))/(A.fromIntegral $ eDepth * eHeight)
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft1D_3r mode) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft1D_r mode) $
#endif
                  fft sign (Z:.depth :.height) width
    in
    case mode of
      Inverse -> A.map (/scale) (go arr)
      _       -> go arr

-- Cube Transform
-- --------------

-- | Discrete Fourier Transform of a 3D array.
--
-- The default implementation requires the array dimensions to be powers of two
-- (else error).
--
fft3D :: FFTElt e
      => Mode
      -> Array DIM3 (Complex e)
      -> Acc (Array DIM3 (Complex e))
fft3D mode arr
  = fft3D' mode (arrayShape arr) (use arr)


-- | Discrete Fourier Transform of a 3D array.
--
-- The default implementation requires the array dimensions to be powers of two.
-- The FFI-backed implementations ignore the Haskell-side size parameter (second
-- argument).
--
fft3D' :: forall e. FFTElt e
       => Mode
       -> DIM3
       -> Acc (Array DIM3 (Complex e))
       -> Acc (Array DIM3 (Complex e))
fft3D' mode (Z :. depth :. height :. width) arr
  = let sign    = signOfMode mode :: e
        scale   = A.fromIntegral (A.size arr)
        go      =
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
                  foreignAcc (Native.fft3D mode) $
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
                  foreignAcc (PTX.fft3D mode) $
#endif
                  fft'

        fft' a  = rotate3D . fft sign (Z:.depth :.height) width
              >-> rotate3D . fft sign (Z:.height:.width)  depth
              >-> rotate3D . fft sign (Z:.width :.depth)  height
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

