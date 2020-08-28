{-# LANGUAGE RankNTypes #-}
-- |
-- Module      : Test.Base
-- Copyright   : [2017..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Test.Base
  where

import Data.Array.Accelerate                                        ( Z(..), (:.)(..), DIM1, DIM2, DIM3, Shape, Elt, Acc, Array )
import Data.Array.Accelerate.Sugar.Array                            ( fromList )
import Data.Array.Accelerate.Sugar.Shape                            ( size )
import Data.Array.Accelerate.Trafo                                  ( Afunction )
import Data.Array.Accelerate.Trafo.Sharing                          ( AfunctionR )
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Math.FFT

import Hedgehog
import qualified Hedgehog.Gen                                       as Gen
import qualified Hedgehog.Range                                     as Range

import Prelude                                                      as P


type RunN = forall f. Afunction f => f -> AfunctionR f

type Transform sh e = Mode -> Acc (Array sh e) -> Acc (Array sh e)


f32 :: Gen Float
f32 = Gen.realFloat (Range.linearFracFrom 0 (-1) 1)

f64 :: Gen Double
f64 = Gen.realFloat (Range.linearFracFrom 0 (-1) 1)

complex :: Gen a -> Gen (Complex a)
complex f = (:+) <$> f <*> f

dim1 :: Gen DIM1
dim1 = (Z :.) <$> Gen.int (Range.linear 1 1024)

dim2 :: Gen DIM2
dim2 = do
  x <- Gen.int (Range.linear 1 128)
  y <- Gen.int (Range.linear 1 48)
  return (Z :. y :. x)

dim3 :: Gen DIM3
dim3 = do
  x <- Gen.int (Range.linear 1 64)
  y <- Gen.int (Range.linear 1 32)
  z <- Gen.int (Range.linear 1 16)
  return (Z :. z :. y :. x)

array :: (Shape sh, Elt e) => sh -> Gen e -> Gen (Array sh e)
array sh gen = fromList sh <$> Gen.list (Range.singleton (size sh)) gen

