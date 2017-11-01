{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE RankNTypes          #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Test.FFT
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Test.FFT ( testFFT )
  where

import Test.ShowType

import Data.Array.Accelerate                                        ( Shape, Slice, Elt, Array, Acc, Exp, (:.)(..), DIM2 )
import Data.Array.Accelerate.Trafo                                  ( Afunction, AfunctionR )
import Data.Array.Accelerate.Array.Sugar                            ( rank )
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT
import Data.Array.Accelerate.Hedgehog.Similar
import qualified Data.Array.Accelerate                              as A
import qualified Data.Array.Accelerate.Hedgehog.Gen.Array           as Gen
import qualified Data.Array.Accelerate.Hedgehog.Gen.Shape           as Gen

import Hedgehog
import qualified Hedgehog.Gen                                       as Gen
import qualified Hedgehog.Range                                     as Range

import Data.Proxy
import Prelude                                                      as P hiding ( reverse )

import Test.Tasty
import Test.Tasty.Hedgehog


-- type Run  = forall a. Arrays a => Acc a -> a
type RunN = forall f. Afunction f => f -> AfunctionR f

type Transform sh e = Mode -> Acc (Array sh e) -> Acc (Array sh e)


floating :: P.RealFloat a => Gen a
floating = Gen.realFloat (Range.linearFracFrom 0 (-1) 1)

complex :: Gen a -> Gen (Complex a)
complex f = (:+) <$> f <*> f

shape :: forall sh. (Gen.Shape sh, Shape sh) => Gen sh
shape = Gen.shape (Range.linear 1 (512 `quot` (2 ^ r)))
  where
    r = rank (undefined::sh)

mode :: Gen Mode
mode = Gen.element [Forward, Reverse, Inverse]

reverse
    :: forall sh e. (Shape sh, Slice sh, Elt e)
    => Acc (Array (sh:.Int) e)
    -> Acc (Array (sh:.Int) e)
reverse arr =
  let sh = A.shape arr
      n  = A.indexHead sh
  in
  A.backpermute sh (\(A.unlift -> ix:.k :: Exp sh :. Exp Int) -> A.lift (ix :. (-k) `mod` n)) arr

norm2
    :: (A.Floating e, Shape sh)
    => Acc (Array (sh:.Int) (Complex e))
    -> Acc (Array sh e)
norm2 = A.map sqrt . A.sum . A.map (\c -> real c * real c + imag c * imag c)

dotc :: (A.RealFloat e, Shape sh)
     => Acc (Array (sh:.Int) (Complex e))
     -> Acc (Array (sh:.Int) (Complex e))
     -> Acc (Array sh (Complex e))
dotc xs ys = A.sum $ A.zipWith (*) xs (A.map conjugate ys)


test_homogeneity
    :: (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform sh (Complex e)
    -> Property
test_homogeneity _ runN transform
  = property
  $ do
        sign  <- forAll mode
        sh    <- forAll shape
        arr   <- forAll (Gen.array sh (complex floating))
        x     <- A.the . A.unit . A.constant <$> forAll (complex floating)
        --
        runN (transform sign . A.map (x *)) arr ~~~ runN (A.map (x *) . transform sign) arr

test_additivity
    :: (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform sh (Complex e)
    -> Property
test_additivity _ runN transform
  = property
  $ do
      sign <- forAll mode
      sh   <- forAll shape
      xs   <- forAll (Gen.array sh (complex floating))
      ys   <- forAll (Gen.array sh (complex floating))
      --
      runN (\u v -> transform sign (A.zipWith (+) u v)) xs ys ~~~ runN (\u v -> A.zipWith (+) (transform sign u) (transform sign v)) xs ys

test_inverse
    :: (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform sh (Complex e)
    -> Property
test_inverse _ runN transform
  = property
  $ do
      sh <- forAll shape
      xs <- forAll (Gen.array sh (complex floating))
      --
      xs ~~~ runN (transform Inverse . transform Forward) xs

test_reverse
    :: (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Slice sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform (sh:.Int) (Complex e)
    -> Property
test_reverse _ runN transform
  = property
  $ do
      sign <- forAll mode
      sh   <- forAll shape
      xs   <- forAll (Gen.array sh (complex floating))
      --
      runN (reverse . transform sign) xs ~~~ runN (transform sign . reverse) xs

test_conjugate
    :: (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Slice sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform (sh:.Int) (Complex e)
    -> Property
test_conjugate _ runN transform
  = property
  $ do
      sign <- forAll mode
      sh   <- forAll shape
      xs   <- forAll (Gen.array sh (complex floating))
      --
      runN (A.map conjugate . transform sign) xs ~~~ runN (transform sign . A.map conjugate . reverse) xs

test_isometry
    :: forall sh e. (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Slice sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform (sh:.Int) (Complex e)
    -> Property
test_isometry _ runN transform
  = property
  $ do
      sign      <- forAll (Gen.element [Forward, Reverse])
      sh@(_:.n) <- forAll shape
      xs        <- forAll (Gen.array sh (complex floating))
      let n'     = A.the . A.unit . A.constant $ sqrt (fromIntegral n)
      --
      runN (norm2 . transform sign) xs ~~~ runN (A.map (n' *) . norm2) xs

test_unitarity
    :: forall sh e. (FFTElt e, Similar e, RealFloat e, Gen.Shape sh, Shape sh, Slice sh, Eq sh)
    => Proxy e
    -> RunN
    -> Transform (sh:.Int) (Complex e)
    -> Property
test_unitarity _ runN transform
  = property
  $ do
      sign      <- forAll (Gen.element [Forward, Reverse])
      sh@(_:.n) <- forAll shape
      xs        <- forAll (Gen.array sh (complex floating))
      ys        <- forAll (Gen.array sh (complex floating))
      let n'     = A.the . A.unit . A.constant $ fromIntegral n
      --
      runN (\u v -> dotc (transform sign u) (transform sign v)) xs ys ~~~ runN (\u v -> A.map (n'*) (dotc u v)) xs ys


testFFT :: RunN -> TestTree
testFFT runN =
  testGroup "FFT"
    [ testFFT' (Proxy::Proxy Float)  runN
    , testFFT' (Proxy::Proxy Double) runN
    ]

testFFT'
    :: forall e. (FFTElt e, Similar e, RealFloat e, Show (ArgType e))
    => Proxy e
    -> RunN
    -> TestTree
testFFT' e runN =
  testGroup (showType e)
    [ testGroup "DIM1"
      [ testProperty "homogeneity" $ test_homogeneity e runN fft1D
      , testProperty "additivity"  $ test_additivity  e runN fft1D
      , testProperty "inverse"     $ test_inverse     e runN fft1D
      , testProperty "reverse"     $ test_reverse     e runN fft1D
      , testProperty "conjugate"   $ test_conjugate   e runN fft1D
      , testProperty "isometry"    $ test_isometry    e runN fft1D
      , testProperty "unitarity"   $ test_unitarity   e runN fft1D
      ]
    , testGroup "DIM2"
      [ testProperty "homogeneity" $ test_homogeneity e runN fft2D
      , testProperty "additivity"  $ test_additivity  e runN fft2D
      , testProperty "inverse"     $ test_inverse     e runN fft2D
      , testProperty "reverse"     $ test_reverse     e runN (fft :: Transform DIM2 (Complex e))
      , testProperty "conjugate"   $ test_conjugate   e runN (fft :: Transform DIM2 (Complex e))
      , testProperty "isometry"    $ test_isometry    e runN (fft :: Transform DIM2 (Complex e))
      , testProperty "unitarity"   $ test_unitarity   e runN (fft :: Transform DIM2 (Complex e))
      ]
    -- , testGroup "DIM3"
    --   [ testProperty "homogeneity" $ test_homogeneity e runN fft3D
    --   , testProperty "additivity"  $ test_additivity  e runN fft3D
    --   , testProperty "inverse"     $ test_inverse     e runN fft3D
    --   , testProperty "reverse"     $ test_reverse     e runN (fft :: Transform DIM3 (Complex e))
    --   , testProperty "conjugate"   $ test_conjugate   e runN (fft :: Transform DIM3 (Complex e))
    --   , testProperty "isometry"    $ test_isometry    e runN (fft :: Transform DIM3 (Complex e))
    --   , testProperty "unitarity"   $ test_unitarity   e runN (fft :: Transform DIM3 (Complex e))
    --   ]
    ]

