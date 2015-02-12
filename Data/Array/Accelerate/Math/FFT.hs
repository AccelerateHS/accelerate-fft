{-# LANGUAGE CPP                      #-}
{-# LANGUAGE EmptyDataDecls           #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE GADTs                    #-}
{-# LANGUAGE ScopedTypeVariables      #-}
{-# LANGUAGE TypeFamilies             #-}
{-# LANGUAGE TypeOperators            #-}
{-# LANGUAGE ViewPatterns             #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT
-- Copyright   : [2012..2013] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
--               [2013..2014] Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Manuel M T Chakravarty <chak@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Computation of a Discrete Fourier Transform using the Cooley-Tuckey
-- algorithm. The time complexity is O(n log n) in the size of the input.
--
-- This uses a naïve divide-and-conquer algorithm whose absolute performance is
-- appalling.
--
module Data.Array.Accelerate.Math.FFT (

  Mode(..),
  fft1D, fft1D',
  fft2D, fft2D',
  fft3D, fft3D',
  fft

) where

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Array.Sugar        ( showShape )
import Data.Array.Accelerate.Data.Complex

#ifdef ACCELERATE_CUDA_BACKEND
import Data.Array.Accelerate.CUDA.Foreign
import Data.Array.Accelerate.Array.Sugar        as S ( shapeToList, shape, EltRepr )
import Data.Array.Accelerate.Type

import Data.Functor
import System.Mem.Weak
import System.IO.Unsafe
import Foreign.CUDA.FFT
import qualified Foreign.CUDA.Types             as CUDA
import qualified Foreign.CUDA.Driver            as CUDA
#endif

import Data.Bits

data Mode = Forward | Reverse | Inverse
  deriving (Eq, Show)

isPow2 :: Int -> Bool
isPow2 x = x .&. (x-1) == 0

signOfMode :: Num a => Mode -> a
signOfMode m
  = case m of
      Forward   -> -1
      Reverse   ->  1
      Inverse   ->  1


-- Vector Transform
-- ----------------
--
-- Discrete Fourier Transform of a vector. Array dimensions must be powers of
-- two else error.
--
fft1D :: (Elt e, IsFloating e)
      => Mode
      -> Vector (Complex e)
      -> Acc (Vector (Complex e))
fft1D mode vec
  = let Z :. len = arrayShape vec
    in
    fft1D' mode len (use vec)

fft1D' :: forall e. (Elt e, IsFloating e)
       => Mode
       -> Int
       -> Acc (Vector (Complex e))
       -> Acc (Vector (Complex e))
fft1D' mode len vec
  = let sign    = signOfMode mode :: e
        scale   = P.fromIntegral len
#ifdef ACCELERATE_CUDA_BACKEND
        sh      = (Z:.len)
        vec'    = cudaFFT mode sh fft' vec
#else
        vec'    = fft' vec
#endif
        fft' a  = fft sign Z len a
    in
    if P.not (isPow2 len)
       then error $ unlines
              [ "Data.Array.Accelerate.FFT: fft1D"
              , "  Array dimensions must be powers of two, but are: " P.++ showShape (Z:.len) ]

       else case mode of
                 Inverse -> A.map (/scale) vec'
                 _       -> vec'


-- Matrix Transform
-- ----------------
--
-- Discrete Fourier Transform of a matrix. Array dimensions must be powers of
-- two else error.
--
fft2D :: (Elt e, IsFloating e)
      => Mode
      -> Array DIM2 (Complex e)
      -> Acc (Array DIM2 (Complex e))
fft2D mode arr
  = let Z :. height :. width = arrayShape arr
    in
    fft2D' mode width height (use arr)


fft2D' :: forall e. (Elt e, IsFloating e)
       => Mode
       -> Int   -- ^ width
       -> Int   -- ^ height
       -> Acc (Array DIM2 (Complex e))
       -> Acc (Array DIM2 (Complex e))
fft2D' mode width height arr
  = let sign    = signOfMode mode :: e
        scale   = P.fromIntegral (width * height)
#ifdef ACCELERATE_CUDA_BACKEND
        sh      = (Z:.height:.width)
        arr'    = cudaFFT mode sh fft' arr
#else
        arr'    = fft' arr
#endif
        fft' a  = A.transpose . fft sign (Z:.width)  height
              >-> A.transpose . fft sign (Z:.height) width
                $ a
    in
    if P.not (isPow2 width && isPow2 height)
       then error $ unlines
              [ "Data.Array.Accelerate.FFT: fft2D"
              , "  Array dimensions must be powers of two, but are: " P.++ showShape (Z:.height:.width) ]

       else case mode of
                 Inverse -> A.map (/scale) arr'
                 _       -> arr'


-- Cube Transform
-- --------------
--
-- Discrete Fourier Transform of a 3D array. Array dimensions must be power of
-- two else error.
--
fft3D :: (Elt e, IsFloating e)
      => Mode
      -> Array DIM3 (Complex e)
      -> Acc (Array DIM3 (Complex e))
fft3D mode arr
  = let Z :. depth :. height :. width = arrayShape arr
    in
    fft3D' mode width height depth (use arr)


fft3D' :: forall e. (Elt e, IsFloating e)
       => Mode
       -> Int   -- ^ width
       -> Int   -- ^ height
       -> Int   -- ^ depth
       -> Acc (Array DIM3 (Complex e))
       -> Acc (Array DIM3 (Complex e))
fft3D' mode width height depth arr
  = let sign    = signOfMode mode :: e
        scale   = P.fromIntegral (width * height)
#ifdef ACCELERATE_CUDA_BACKEND
        sh      = (Z:.depth:.height:.width)
        arr'    = cudaFFT mode sh fft' arr
#else
        arr'    = fft' arr
#endif
        fft' a  = rotate3D . fft sign (Z:.width :.depth)  height
              >-> rotate3D . fft sign (Z:.height:.width)  depth
              >-> rotate3D . fft sign (Z:.depth :.height) width
                $ a
    in
    if P.not (isPow2 width && isPow2 height && isPow2 depth)
       then error $ unlines
              [ "Data.Array.Accelerate.FFT: fft3D"
              , "  Array dimensions must be powers of two, but are: " P.++ showShape (Z:.depth:.height:.width) ]

       else case mode of
                 Inverse -> A.map (/scale) arr'
                 _       -> arr'



rotate3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
rotate3D arr
  = backpermute (swap (A.shape arr)) swap arr
  where
    swap :: Exp DIM3 -> Exp DIM3
    swap ix =
      let Z :. m :. k :. l = unlift ix  :: Z :. Exp Int :. Exp Int :. Exp Int
      in  lift $ Z :. k :. l :. m


-- Rank-generalised Cooley-Tuckey DFT
--
-- We require the innermost dimension be passed as a Haskell value because we
-- can't do divide-and-conquer recursion directly in the meta-language.
--
fft :: forall sh e. (Slice sh, Shape sh, IsFloating e, Elt e)
    => e
    -> sh
    -> Int
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fft sign sh sz arr = go sz 0 1
  where
    go :: Int -> Int -> Int -> Acc (Array (sh:.Int) (Complex e))
    go len offset stride
      | len == 2
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
          sz' ==* 0 ? ( (arr ! lift (sh' :. offset')) + (arr ! lift (sh' :. offset' + stride'))
          {-  ==* 1-} , (arr ! lift (sh' :. offset')) - (arr ! lift (sh' :. offset' + stride')) )

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

#ifdef ACCELERATE_CUDA_BACKEND
-- FFT using the CUFFT library to enable high performance for the CUDA backend of
-- Accelerate. The implementation works on all arrays of rank less than or equal
-- to 3. The result is un-normalised.
--
cudaFFT :: forall e sh. (Shape sh, Elt e, IsFloating e)
        => Mode
        -> sh
        -> (Acc (Array sh (Complex e)) -> Acc (Array sh (Complex e)))
        -> Acc (Array sh (Complex e))
        -> Acc (Array sh (Complex e))
cudaFFT mode sh = cudaFFT'
  where
    -- Plan the FFT.
    -- Doing this in unsafePerformIO so it is not reperformed every time the
    -- AST is evaluated.
    --
    hndl = unsafePerformIO $ do
            plan <- case shapeToList sh of
                     [width]                -> plan1D              width types 1
                     [width, height]        -> plan2D       height width types
                     [width, height, depth] -> plan3D depth height width types
                     _                      -> error "Accelerate-fft cannot use CUFFT for arrays of dimensions higher than 3"
            addFinalizer plan (destroy plan)
            return plan

    types = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> C2C
              TypeDouble{}  -> Z2Z
              TypeCFloat{}  -> C2C
              TypeCDouble{} -> Z2Z

    cudaFFT' p arr = deinterleave sh (foreignAcc ff pure (interleave arr))
      where
        ff          = CUDAForeignAcc "foreignFFT" foreignFFT
        -- Unfortunately the pure version of the function needs to be wrapped in
        -- interleave and deinterleave to match how the foreign version works.
        --
        -- RCE: Do the interleaving and deinterleaving in foreignFFT
        --
        -- TLM: The interleaving might get fused into other parts of the
        --      computation and thus be okay. We should really support multi types
        --      such as float2 instead.
        --
        pure        = interleave . p . deinterleave sh
        sign        = signOfMode mode :: Int

        foreignFFT :: CUDA.Stream -> Array DIM1 e -> CIO (Array DIM1 e)
        foreignFFT stream arr' = do
          output <- allocateArray (S.shape arr')
          iptr   <- floatingDevicePtr arr'
          optr   <- floatingDevicePtr output

          --Execute
          liftIO $ do
            setStream hndl stream
            execute iptr optr

          return output

        execute :: CUDA.DevicePtr e -> CUDA.DevicePtr e -> IO ()
        execute iptr optr
          = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> execC2C hndl iptr optr sign
              TypeDouble{}  -> execZ2Z hndl iptr optr sign
              TypeCFloat{}  -> execC2C hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign
              TypeCDouble{} -> execZ2Z hndl (CUDA.castDevPtr iptr) (CUDA.castDevPtr optr) sign

        floatingDevicePtr :: Vector e -> CIO (CUDA.DevicePtr e)
        floatingDevicePtr v
          = case (floatingType :: FloatingType e) of
              TypeFloat{}   -> singleDevicePtr v
              TypeDouble{}  -> singleDevicePtr v
              TypeCFloat{}  -> CUDA.castDevPtr <$> singleDevicePtr v
              TypeCDouble{} -> CUDA.castDevPtr <$> singleDevicePtr v

        singleDevicePtr :: DevicePtrs (EltRepr e) ~ CUDA.DevicePtr b => Vector e -> CIO (CUDA.DevicePtr b)
        singleDevicePtr v = devicePtrsOfArray v
#endif

-- Append two arrays. Doesn't do proper bounds checking or intersection...
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
                     in  i <* n ? (xs ! lift (sz:.i), ys ! lift (sz:.i-n) ))


#ifdef ACCELERATE_CUDA_BACKEND
{-# RULES
  "interleave/deinterleave" forall sh x. deinterleave sh (interleave x) = x;
  "deinterleave/interleave" forall sh x. interleave (deinterleave sh x) = x
 #-}

-- Interleave the real and imaginary components in a complex array and produce a
-- flattened vector. This allows us to mimic the float2 structure used by CUFFT
-- to store complex numbers.
--
{-# NOINLINE interleave #-}
interleave :: (Shape sh, Elt e) => Acc (Array sh (Complex e)) -> Acc (Vector e)
interleave arr = generate sh swizzle
  where
    sh          = index1 (2 * A.size arr)
    swizzle ix  =
      let i = indexHead ix
          v = arr A.!! (i `div` 2)
      in
      i `mod` 2 ==* 0 ? (real v, imag v)

-- Deinterleave a vector into a complex array. Assumes the array is even in length.
--
{-# NOINLINE deinterleave #-}
deinterleave :: (Shape sh, Elt e) => sh -> Acc (Vector e) -> Acc (Array sh (Complex e))
deinterleave (constant -> sh) arr =
  generate sh (\ix -> let i = toIndex sh ix * 2
                      in  lift (arr A.!! i :+ arr A.!! (i+1)))
#endif
