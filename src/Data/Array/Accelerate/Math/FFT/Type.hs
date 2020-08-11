{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs             #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# OPTIONS_HADDOCK hide #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Type
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.Type
  where

import Data.Array.Accelerate                                        as A
import Data.Array.Accelerate.Data.Complex                           as A
import Data.Array.Accelerate.Sugar.Elt
import Data.Primitive.Vec


-- For explicit dictionary reification, to discover the concrete type the
-- operation should be performed at.
--
data NumericR a where
  NumericRfloat32 :: NumericR Float
  NumericRfloat64 :: NumericR Double

class (RealFloat a, FromIntegral Int a, Elt a, EltR (Complex a) ~ Vec2 a) => Numeric a where
  numericR :: NumericR a

instance Numeric Float where
  numericR = NumericRfloat32

instance Numeric Double where
  numericR = NumericRfloat64

