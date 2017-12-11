{-# LANGUAGE BangPatterns        #-}
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

import Test.Base
import Test.ShowType

import Data.Array.Accelerate                                        as A hiding ( RealFloat, Eq, reverse )
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT
import Data.Array.Accelerate.Test.Similar
import qualified Data.Array.Accelerate                              as A

import Hedgehog
import qualified Hedgehog.Gen                                       as Gen

import Test.Tasty
import Test.Tasty.Hedgehog

import Data.Proxy
import Prelude                                                      as P hiding ( reverse )


testFFT :: RunN -> TestTree
testFFT runN =
  testGroup "FFT"
    [ testFFT' f32 runN
    , testFFT' f64 runN
    ]

testFFT'
    :: forall e. (FFTElt e, Similar e, RealFloat e, Show (ArgType e))
    => Gen e
    -> RunN
    -> TestTree
testFFT' e runN =
  testGroup (showType (Proxy::Proxy e))
    [ testGroup "DIM1"
      [ testProperty "homogeneity" $ test_homogeneity runN fft1D dim1 e
      , testProperty "additivity"  $ test_additivity  runN fft1D dim1 e
      , testProperty "inverse"     $ test_inverse     runN fft1D dim1 e
      , testProperty "reverse"     $ test_reverse     runN fft1D dim1 e
      , testProperty "conjugate"   $ test_conjugate   runN fft1D dim1 e
      , testProperty "isometry"    $ test_isometry    runN fft1D dim1 e
      , testProperty "unitarity"   $ test_unitarity   runN fft1D dim1 e
      ]
    , testGroup "DIM2"
      [ testProperty "homogeneity" $ test_homogeneity runN fft2D dim2 e
      , testProperty "additivity"  $ test_additivity  runN fft2D dim2 e
      , testProperty "inverse"     $ test_inverse     runN fft2D dim2 e
      , testProperty "reverse"     $ test_reverse     runN fft   dim2 e
      , testProperty "conjugate"   $ test_conjugate   runN fft   dim2 e
      , testProperty "isometry"    $ test_isometry    runN fft   dim2 e
      , testProperty "unitarity"   $ test_unitarity   runN fft   dim2 e
      ]
    , testGroup "DIM3"
      [ testProperty "homogeneity" $ test_homogeneity runN fft3D dim3 e
      , testProperty "additivity"  $ test_additivity  runN fft3D dim3 e
      , testProperty "inverse"     $ test_inverse     runN fft3D dim3 e
      , testProperty "reverse"     $ test_reverse     runN fft   dim3 e
      , testProperty "conjugate"   $ test_conjugate   runN fft   dim3 e
      , testProperty "isometry"    $ test_isometry    runN fft   dim3 e
      , testProperty "unitarity"   $ test_unitarity   runN fft   dim3 e
      ]
    ]


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

scalar :: Elt e => e -> Scalar e
scalar x = fromFunction Z (const x)


test_homogeneity
    :: (FFTElt e, Similar e, Shape sh, Eq sh)
    => RunN
    -> Transform sh (Complex e)
    -> Gen sh
    -> Gen e
    -> Property
test_homogeneity runN transform dim e =
  property $ do
    sign  <- forAll mode
    sh    <- forAll dim
    arr   <- forAll (array sh (complex e))
    x     <- forAll (complex e)
    --
    let !go1 = runN (\u -> transform sign . A.map (the u *))
        !go2 = runN (\u -> A.map (the u *) . transform sign)
    --
    go1 (scalar x) arr ~~~ go2 (scalar x) arr

test_additivity
    :: (FFTElt e, Similar e, Shape sh, Eq sh)
    => RunN
    -> Transform sh (Complex e)
    -> Gen sh
    -> Gen e
    -> Property
test_additivity runN transform dim e =
  property $ do
    sign <- forAll mode
    sh   <- forAll dim
    xs   <- forAll (array sh (complex e))
    ys   <- forAll (array sh (complex e))
    --
    let !go1 = runN (\u v -> transform sign (A.zipWith (+) u v))
        !go2 = runN (\u v -> A.zipWith (+) (transform sign u) (transform sign v))
    --
    go1 xs ys ~~~ go2 xs ys

test_inverse
    :: (FFTElt e, Similar e, Shape sh, Eq sh)
    => RunN
    -> Transform sh (Complex e)
    -> Gen sh
    -> Gen e
    -> Property
test_inverse runN transform dim e =
  property $ do
    sh <- forAll dim
    xs <- forAll (array sh (complex e))
    --
    let !go = runN (transform Inverse . transform Forward)
    xs ~~~ go xs

test_reverse
    :: (FFTElt e, Similar e, Shape sh, Slice sh, Eq sh)
    => RunN
    -> Transform (sh:.Int) (Complex e)
    -> Gen (sh:.Int)
    -> Gen e
    -> Property
test_reverse runN transform dim e =
  property $ do
    sign <- forAll mode
    sh   <- forAll dim
    xs   <- forAll (array sh (complex e))
    --
    let !go1 = runN (reverse . transform sign)
        !go2 = runN (transform sign . reverse)
    --
    go1 xs ~~~ go2 xs

test_conjugate
    :: (FFTElt e, Similar e, Shape sh, Slice sh, Eq sh)
    => RunN
    -> Transform (sh:.Int) (Complex e)
    -> Gen (sh:.Int)
    -> Gen e
    -> Property
test_conjugate runN transform dim e =
  property $ do
    sign <- forAll mode
    sh   <- forAll dim
    xs   <- forAll (array sh (complex e))
    --
    let !go1 = runN (A.map conjugate . transform sign)
        !go2 = runN (transform sign . A.map conjugate . reverse)
    --
    go1 xs ~~~ go2 xs

test_isometry
    :: forall sh e. (FFTElt e, Similar e, Shape sh, Slice sh, Eq sh)
    => RunN
    -> Transform (sh:.Int) (Complex e)
    -> Gen (sh:.Int)
    -> Gen e
    -> Property
test_isometry runN transform dim e =
  property $ do
    sign      <- forAll (Gen.element [Forward, Reverse])
    sh@(_:.n) <- forAll dim
    xs        <- forAll (array sh (complex e))
    --
    let !go1   = runN (norm2 . transform sign)
        !go2   = runN (\u -> A.map (the u *) . norm2)
    --
    go1 xs ~~~ go2 (scalar (sqrt (P.fromIntegral n))) xs

test_unitarity
    :: forall sh e. (FFTElt e, Similar e, RealFloat e, Shape sh, Slice sh, Eq sh)
    => RunN
    -> Transform (sh:.Int) (Complex e)
    -> Gen (sh:.Int)
    -> Gen e
    -> Property
test_unitarity runN transform dim e =
  property $ do
    sign      <- forAll (Gen.element [Forward, Reverse])
    sh@(_:.n) <- forAll dim
    xs        <- forAll (array sh (complex e))
    ys        <- forAll (array sh (complex e))
    --
    let !go1   = runN (\u v -> dotc (transform sign u) (transform sign v))
        !go2   = runN (\m u v -> A.map (the m *) (dotc u v))
    --
    go1 xs ys ~~~ go2 (scalar (P.fromIntegral n)) xs ys

