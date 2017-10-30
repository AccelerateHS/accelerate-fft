{-# LANGUAGE CPP                 #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE RebindableSyntax    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Adhoc
-- Copyright   : [2017] Henning Thielemann
--               [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Implementation of ad-hoc FFT stolen from the accelerate-fourier by Henning
-- Thielemann (BSD3 licensed), and updated to work with current Accelerate. That
-- package contains other more sophisticated algorithms as well.
--

module Data.Array.Accelerate.Math.FFT.Adhoc ( fft )
  where

import Data.Array.Accelerate                                        hiding ( transpose )
import Data.Array.Accelerate.Data.Bits
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT.Mode
import Data.Array.Accelerate.Math.FFT.Elt

#if __GLASGOW_HASKELL__ < 800
import Prelude                                                      ( fromInteger )
#endif


fft :: (Shape sh, Slice sh, FFTElt e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
fft mode arr =
  let len             = indexHead (shape arr)
      (pow2, smooth5) = is2or5smooth len
  in
  if len <= 1 then arr                        else
  if pow2     then ditSplitRadixLoop mode arr else
  if smooth5  then dit235            mode arr
              else transformChirp235 mode arr


-- Implementations
-- ---------------

is2or5smooth :: Exp Int -> (Exp Bool, Exp Bool)
is2or5smooth len =
  let maxPowerOfTwo = len .&. negate len
      lenOdd        = len `quot` maxPowerOfTwo
  in
  ( 1 == lenOdd
  , 1 == divideMaxPower 5 (divideMaxPower 3 lenOdd)
  )

divideMaxPower :: Exp Int -> Exp Int -> Exp Int
divideMaxPower fac =
  while (\n -> n `rem`  fac == 0)
        (\n -> n `quot` fac)


-- -- | Split-radix for power-of-two sizes
-- --
-- ditSplitRadix
--     :: (Shape sh, Slice sh, FFTElt e)
--     => Mode
--     -> Acc (Array (sh:.Int) (Complex e))
--     -> Acc (Array (sh:.Int) (Complex e))
-- ditSplitRadix mode arr =
--   if indexHead (shape arr) <= 1
--     then arr
--     else ditSplitRadixLoop mode arr

ditSplitRadixLoop
    :: forall sh e. (Shape sh, Slice sh, FFTElt e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
ditSplitRadixLoop mode arr =
  let
      twiddleSR (fromIntegral -> n4) k (fromIntegral -> j) =
        let w = pi * k * j / (2*n4)
        in  lift (cos w :+ signOfMode mode * sin w)

      twiddle len4 k =
        generate (index1 len4) (twiddleSR len4 k . indexHead)

      step (unlift -> (us,zs)) =
        let
            k           = indexHead (shape zs)
            tw1         = twiddle k 1
            tw3         = twiddle k 3
            --
            im          = constant (0 :+ signOfMode mode)
            twidZeven   = zipWithExtrude1 (*) tw1 (sieveV 2 0 zs)
            twidZodd    = zipWithExtrude1 (*) tw3 (sieveV 2 1 zs)
            zsum        = zipWith (+) twidZeven twidZodd
            zdiff       = map (im *) (zipWith (-) twidZeven twidZodd)
            zcomplete   = zsum ++ zdiff
            _ :. n :. _ = unlift (shape zcomplete) :: Exp sh :. Exp Int :. Exp Int
        in
        lift ( zipWith (+) us zcomplete ++ zipWith (-) us zcomplete
             , dropV n us
             )

      rebase s = lift (transform2 (-1) (afst s), asnd s)

      reorder (unlift -> (xs,ys)) =
        let evens = sieve 2 0 xs
            odds  = sieve 2 1 xs
        in
        lift (evens ++^ ys, twist 2 odds)

      initial =
        let sh :. n = unlift (shape arr) :: Exp sh :. Exp Int
        in  lift ( reshape (lift (sh :. constant 1 :. n)) arr
                 , fill    (lift (sh :. constant 0 :. n `quot` 2)) 0
                 )
  in
  headV
    $ afst
    $ awhile (\s -> unit (indexHead (indexTail (shape (asnd s))) > 0)) step
    $ rebase
    $ awhile (\s -> unit (indexHead (shape (asnd s)) > 1)) reorder
    $ initial


-- | Decimation in time for sizes that are composites of the factors 2,3 and 5.
-- These sizes are known as 5-smooth numbers or the Hamming sequence.
--
-- <http://oeis.org/A051037>
--
dit235
    :: forall sh e. (Shape sh, Slice sh, FFTElt e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
dit235 mode arr =
  let
      merge :: forall sh' a. (Shape sh', Slice sh', Elt a)
            => Acc (Array (sh':.Int:.Int) a)
            -> Acc (Array (sh':.Int) a)
      merge xs =
        let sh :. m :. n = unlift (shape xs) :: Exp sh' :. Exp Int :. Exp Int
        in  backpermute
              (lift (sh :. m*n))
              (\(unlift -> ix :. k :: Exp sh' :. Exp Int) ->
                  let (q,r) = k `quotRem` m
                  in  lift (ix :. r :. q))
              xs

      step fac xs =
        let sh :. count :. len = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
            twiddled           = transpose
                               $ zipWithExtrude2 (*) (twiddleFactors fac len)
                               $ reshape (lift (sh :. count `quot` fac :. fac :. len)) xs
        in
        merge $ if fac == 5 then transform5 cache5 twiddled else
                if fac == 4 then transform4 cache4 twiddled else
                if fac == 3 then transform3 cache3 twiddled
                            else transform2 cache2 twiddled

      initial :: Acc (Array (sh:.Int:.Int) (Complex e), Vector Int)
      initial =
        let sh :. n = unlift (shape arr) :: Exp sh :. Exp Int
        in  lift ( reshape (lift (sh :. constant 1 :. n)) arr
                 , fill (index1 0) 0
                 )

      twiddleFactors :: Exp Int -> Exp Int -> Acc (Matrix (Complex e))
      twiddleFactors m n =
        generate (index2 m n)
                 (\(unlift -> Z :. j :. i) -> twiddle (m*n) j i)

      cisrat :: Exp Int -> Exp Int -> Exp (Complex e)
      cisrat d n =
        let w = 2*pi * fromIntegral n / fromIntegral d
        in  lift (cos w :+ signOfMode mode * sin w)

      twiddle :: Exp Int -> Exp Int -> Exp Int -> Exp (Complex e)
      twiddle n k j = cisrat n ((k*j) `rem` n)

      cache2 :: Exp (Complex e)
      cache2 = -1

      cache3 :: Exp (Complex e, Complex e)
      cache3 =
        let sqrt3d2 = sqrt 3 / 2
            mhalf   = constant (-1/2)
            s       = signOfMode mode
            u       = s * sqrt3d2
        in
        lift (mhalf :+ u, mhalf :+ (-u))

      cache4 :: Exp (Complex e, Complex e, Complex e)
      cache4 =
        let s = signOfMode mode
        in  lift (constant 0 :+ s, constant ((-1) :+ (-0)), constant 0 :+ (-s))

      cache5 :: Exp (Complex e, Complex e, Complex e, Complex e)
      cache5 =
        let z = cisrat 5
        in  lift (z 1, z 2, z 3, z 4)
  in
  headV
    $ afst
    $ awhile
        (\s -> unit (length (asnd s) > 0))
        (\s -> let (xs,fs) = unlift s
                   f       = fs !! 0
               in
               lift (step f xs, tail fs))
    $ awhile
        (\s -> unit (indexHead (shape (afst s)) > 1))
        (\s -> let (xs,fs)      = unlift s
                   len          = indexHead (shape xs)
                   divides k n  = n `rem` k == 0
                   f            = if divides 3 len then 3 else
                                  if divides 4 len then 4 else
                                  if divides 5 len then 5
                                                   else 2
               in
               lift (twist f xs, unit f `cons` fs))
    $ initial


-- | Transformation of arbitrary length base on Bluestein on a 5-smooth size.
--
transformChirp235
    :: (Shape sh, Slice sh, FFTElt e)
    => Mode
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
transformChirp235 mode arr =
  let n = indexHead (shape arr)
      f = ceiling5Smooth (2*n)
  in
  transformChirp mode f (dit235 Forward) (dit235 Inverse) arr


transformChirp
    :: (Shape sh, Slice sh, FFTElt e)
    => Mode
    -> Exp Int
    -> (forall sh'. (Shape sh', Slice sh') => Acc (Array (sh':.Int) (Complex e)) -> Acc (Array (sh':.Int) (Complex e)))
    -> (forall sh'. (Shape sh', Slice sh') => Acc (Array (sh':.Int) (Complex e)) -> Acc (Array (sh':.Int) (Complex e)))
    -> Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array (sh:.Int) (Complex e))
transformChirp mode p analysis synthesis arr =
  let sz :. n   = unlift (shape arr)
      --
      chirp     =
        generate (index1 p) $ \ix ->
          let k  = unindex1 ix
              sk = fromIntegral (if p > 2*k then k else k-p)
              w  = pi * sk * sk / fromIntegral n
          in
          lift $ cos w :+ signOfMode mode * sin w
      --
      spectrum  = analysis
                $ map conjugate chirp
                  `consV`
                  reshape (lift (Z :. shapeSize sz :. p))
                          (pad p 0 (zipWithExtrude1 (*) chirp arr))
      scaleDown xs =
        let scale x (unlift -> r :+ i) = lift (x*r :+ x*i)
            len                        = indexHead (shape xs)
        in  map (scale (recip (fromIntegral len))) xs
  in
  if n <= 1
    then arr
    else take n
       $ scaleDown
       $ zipWithExtrude1 (*) chirp
       $ synthesis
       $ zipWithExtrude1 (*) (headV spectrum)
       $ reshape (lift (sz:.p)) (tailV spectrum)


ceiling5Smooth :: Exp Int -> Exp Int
ceiling5Smooth n =
  let (i2,i3,i5) = unlift (snd (ceiling5Smooth' (fromIntegral n :: Exp Double)))
  in  pow i2 2 * pow i3 3 * pow i5 5

ceiling5Smooth'
    :: (RealFloat a, Ord a, FromIntegral Int a)
    => Exp a
    -> Exp (a, (Int,Int,Int))
ceiling5Smooth' n =
  let d3 = ceiling (logBase 3 n)
      d5 = ceiling (logBase 5 n)
      --
      argmin x y = if fst x < fst y then x else y
  in
  the $ fold1All argmin
      $ generate (index2 d5 d3) -- this is probably quite small!
                 (\(unlift -> Z :. i5 :. i3) ->
                    let
                        p53 = 5 ** fromIntegral i5 * 3 ** fromIntegral i3
                        i2  = 0 `max` ceiling (logBase 2 (n/p53))
                    in
                    lift ( p53 * 2 ** fromIntegral i2
                         , (i2,i3,i5)
                         ))

-- Utilities
-- ---------

pow :: Exp Int -> Exp Int -> Exp Int
pow x k
  = snd
  $ while (\ip -> fst ip < k)
          (\ip -> lift (fst ip + 1, snd ip * x))
          (constant (0,1))

pad :: (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp e
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
pad n x xs =
  let sz = indexTail (shape xs)
      sh = lift (sz :. n)
  in
  xs ++ fill sh x

cons :: forall sh e. (Shape sh, Slice sh, Elt e)
     => Acc (Array sh e)
     -> Acc (Array (sh:.Int) e)
     -> Acc (Array (sh:.Int) e)
cons x xs =
  let x' = reshape (lift (shape x :. constant 1)) x
  in  x' ++ xs

consV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
consV x xs =
  let sh :. n = unlift (shape x) :: Exp sh :. Exp Int
  in  reshape (lift (sh :. constant 1 :. n)) x ++^ xs

headV :: (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int) e)
headV xs = slice xs (lift (Any :. (0::Int) :. All))

tailV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
tailV = dropV 1

dropV :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Exp Int
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
dropV u xs =
  let sh :. m :. n = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
  in
  backpermute (lift (sh :. 0 `max` (m-u) :. n))
              (\(unlift -> ix :. j :. i :: Exp sh :. Exp Int :. Exp Int) -> lift (ix :. j+u :. i))
              xs

sieve
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp Int
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
sieve fac start xs =
  let sh :. n = unlift (shape xs) :: Exp sh :. Exp Int
  in
  backpermute
    (lift (sh :. n `quot` fac))
    (\(unlift -> ix :. j :: Exp sh :. Exp Int) -> lift (ix :. fac*j + start))
    xs

sieveV
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Exp Int
    -> Exp Int
    -> Acc (Array (sh:.Int:.Int) e)
    -> Acc (Array (sh:.Int:.Int) e)
sieveV fac start xs =
  let sh :. m :. n = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
  in
  backpermute
    (lift (sh :. m `quot` fac :. n))
    (\(unlift -> ix :. j :. i :: Exp sh :. Exp Int :. Exp Int) -> lift (ix :. fac*j+start :. i))
    xs

twist :: forall sh e. (Shape sh, Slice sh, Elt e)
      => Exp Int
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
twist fac xs =
  let sh :. m :. n = unlift (shape xs) :: Exp sh :. Exp Int :. Exp Int
  in
  backpermute
    (lift (sh :. fac*m :. n `quot` fac))
    (\(unlift -> ix :. j :. i :: Exp sh :. Exp Int :. Exp Int) -> lift (ix :. j `quot` fac :. fac*i + j `rem` fac))
    xs


infixr 5 ++^
(++^) :: forall sh e. (Slice sh, Shape sh, Elt e)
      => Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
      -> Acc (Array (sh:.Int:.Int) e)
(++^) xs ys =
  let sh1 :. m1 :. n1 = unlift (shape xs)
      sh2 :. m2 :. n2 = unlift (shape ys)
  in
  generate (lift (intersect sh1 sh2 :. m1+m2 :. min n1 n2))
           (\ix -> let sh :. j :. i = unlift ix :: Exp sh :. Exp Int :. Exp Int
                   in  if j < m1
                         then xs ! ix
                         else ys ! lift (sh :. j-m1 :. i))

zipWithExtrude1
    :: (Shape sh, Slice sh, Elt a, Elt b, Elt c)
    => (Exp a -> Exp b -> Exp c)
    -> Acc (Array DIM1      a)
    -> Acc (Array (sh:.Int) b)
    -> Acc (Array (sh:.Int) c)
zipWithExtrude1 f xs ys =
  zipWith f (replicate (lift (indexTail (shape ys) :. All)) xs) ys

zipWithExtrude2
    :: (Shape sh, Slice sh, Elt a, Elt b, Elt c)
    => (Exp a -> Exp b -> Exp c)
    -> Acc (Array DIM2           a)
    -> Acc (Array (sh:.Int:.Int) b)
    -> Acc (Array (sh:.Int:.Int) c)
zipWithExtrude2 f xs ys =
  zipWith f (replicate (lift (indexTail (indexTail (shape ys)) :. All :. All)) xs) ys

transpose
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Acc (Array (sh:.Int:.Int) e)
    -> Acc (Array (sh:.Int:.Int) e)
transpose arr =
  let swap (unlift -> (ix:.r:.c) :: Exp sh :. Exp Int :. Exp Int) = lift (ix:.c:.r)
  in  backpermute (swap (shape arr)) swap arr

transform2
    :: (Shape sh, Slice sh, Num e)
    => Exp e
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform2 v xs =
  generate
    (lift (indexTail (shape xs) :. (2::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
        in
        if k == 0 then x0+x1
                  else x0+v*x1)

transform3
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform3 (unlift -> (z1,z2)) xs =
  generate
    (lift (indexTail (shape xs) :. (3::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            --
            ((s,_), (zx1,zx2)) = sumAndConvolve2 (x1,x2) (z1,z2)
        in
        if k == 0    then x0 + s   else
        if k == 1    then x0 + zx1
        {- k == 2 -} else x0 + zx2)

transform4
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform4 (unlift -> (z1,z2,z3)) xs =
  generate
    (lift (indexTail (shape xs) :. (4::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            x3 = xs ! lift (ix :. (3::Int))
            --
            x02a = x0+x2
            x02b = x0+z2*x2
            x13a = x1+x3
            x13b = x1+z2*x3
        in
        if k == 0    then x02a +      x13a else
        if k == 1    then x02b + z1 * x13b else
        if k == 2    then x02a + z2 * x13a
        {- k == 3 -} else x02b + z3 * x13b)

-- Use Rader's trick for mapping the transform to a convolution and apply
-- Karatsuba's trick at two levels (i.e. total three times) to that convolution.
--
-- 0 0 0 0 0
-- 0 1 2 3 4
-- 0 2 4 1 3
-- 0 3 1 4 2
-- 0 4 3 2 1
--
-- Permutation.T: 0 1 2 4 3
--
-- 0 0 0 0 0
-- 0 1 2 4 3
-- 0 2 4 3 1
-- 0 4 3 1 2
-- 0 3 1 2 4
--
transform5
    :: forall sh e. (Shape sh, Slice sh, Num e)
    => Exp (e,e,e,e)
    -> Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
transform5 (unlift -> (z1,z2,z3,z4)) xs =
  generate
    (lift (indexTail (shape xs) :. (5::Int)))
    (\(unlift -> ix :. k :: Exp sh :. Exp Int) ->
        let
            x0 = xs ! lift (ix :. (0::Int))
            x1 = xs ! lift (ix :. (1::Int))
            x2 = xs ! lift (ix :. (2::Int))
            x3 = xs ! lift (ix :. (3::Int))
            x4 = xs ! lift (ix :. (4::Int))
            --
            ((s,_), (d1,d2,d4,d3)) = sumAndConvolve4 (x1,x3,x4,x2) (z1,z2,z4,z3)
        in
        if k == 0    then x0 + s  else
        if k == 1    then x0 + d1 else
        if k == 2    then x0 + d2 else
        if k == 3    then x0 + d3
        {- k == 4 -} else x0 + d4)


-- Some small size convolutions using the Karatsuba trick.
--
-- This does not use Toom-3 multiplication, because this requires division by
-- 2 and 6, and thus 'Fractional' constraints.
--
sumAndConvolve2
    :: Num e
    => (Exp e, Exp e)
    -> (Exp e, Exp e)
    -> ((Exp e, Exp e), (Exp e, Exp e))
sumAndConvolve2 (a0,a1) (b0,b1) =
  let sa01   = a0+a1
      sb01   = b0+b1
      ab0ab1 = a0*b0+a1*b1
  in
  ((sa01, sb01), (ab0ab1, sa01*sb01-ab0ab1))

-- sumAndConvolve3
--     :: Num e
--     => (Exp e, Exp e, Exp e)
--     -> (Exp e, Exp e, Exp e)
--     -> ((Exp e, Exp e), (Exp e, Exp e, Exp e))
-- sumAndConvolve3 (a0,a1,a2) (b0,b1,b2) =
--   let ab0   = a0*b0
--       dab12 = a1*b1 - a2*b2
--       sa01  = a0+a1; sb01 = b0+b1; tab01 = sa01*sb01 - ab0
--       sa02  = a0+a2; sb02 = b0+b2; tab02 = sa02*sb02 - ab0
--       sa012 = sa01+a2
--       sb012 = sb01+b2
--       --
--       d0    = sa012*sb012 - tab01 - tab02
--       d1    = tab01 - dab12
--       d2    = tab02 + dab12
--   in
--   ((sa012, sb012), (d0, d1, d2))

sumAndConvolve4
  :: Num e
  => (Exp e, Exp e, Exp e, Exp e)
  -> (Exp e, Exp e, Exp e, Exp e)
  -> ((Exp e, Exp e), (Exp e, Exp e, Exp e, Exp e))
sumAndConvolve4 (a0,a1,a2,a3) (b0,b1,b2,b3) =
  let ab0    = a0*b0
      ab1    = a1*b1
      sa01   = a0+a1; sb01 = b0+b1
      ab01   = sa01*sb01 - (ab0+ab1)
      ab2    = a2*b2
      ab3    = a3*b3
      sa23   = a2+a3; sb23 = b2+b3
      ab23   = sa23*sb23 - (ab2+ab3)
      c0     = ab0  + ab2 - (ab1 + ab3)
      c1     = ab01 + ab23
      ab02   = (a0+a2)*(b0+b2)
      ab13   = (a1+a3)*(b1+b3)
      sa0123 = sa01+sa23
      sb0123 = sb01+sb23
      ab0123 = sa0123*sb0123 - (ab02+ab13)
      --
      d0     = ab13   + c0
      d1     = c1
      d2     = ab02   - c0
      d3     = ab0123 - c1
  in
  ((sa0123, sb0123), (d0, d1, d2, d3))

