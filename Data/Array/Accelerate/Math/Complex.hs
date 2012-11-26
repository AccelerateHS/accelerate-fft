{-# LANGUAGE FlexibleInstances    #-}
{-# LANGUAGE IncoherentInstances  #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# OPTIONS -fno-warn-orphans #-}

module Data.Array.Accelerate.Math.Complex
  where

import Prelude
import Data.Function
import Data.Array.Accelerate                    as A

type Complex a = (a, a)

instance (Elt a, IsFloating a) => Num (Exp (Complex a)) where
  c1 + c2       = lift ( on (+) real c1 c2, on (+) imag c1 c2 )
  c1 - c2       = lift ( on (-) real c1 c2, on (-) imag c1 c2 )
  c1 * c2       = let (x,  y)   = unlift c1
                      (x', y')  = unlift c2     :: Complex (Exp a)
                  in lift (x*x'-y*y', x*y'+y*x')

  negate c      = lift ( negate (real c), negate (imag c) )
  abs z         = lift ( magnitude z, constant 0 )
  signum z      = let r         = magnitude z
                      (x, y)    = unlift z
                  in r ==* 0 ? (constant (0,0), lift (x/r, y/r))

  fromInteger n
    = lift (constant (fromInteger n), constant 0)


instance (Elt a, IsFloating a) => Fractional (Exp (Complex a)) where
  c1 / c2
    = let (a,b) = unlift c1
          (c,d) = unlift c2     :: Complex (Exp a)
          den   = c^(2 :: Int) + d^(2 :: Int)
          re    = (a * c + b * d) / den
          im    = (b * c - a * d) / den
      in
      lift (re, im)

  fromRational x
    = lift (constant (fromRational x), constant 0)


-- | Non-negative magnitude of a complex number
--
magnitude :: (Elt a, IsFloating a) => Exp (Complex a) -> Exp a
magnitude c =
  let (r, i) = unlift c
  in sqrt (r*r + i*i)

-- | The phase of a complex number, in the range (-pi, pi]. If the magnitude is
-- zero, then so is the phase.
--
phase :: (Elt a, IsFloating a) => Exp (Complex a) -> Exp a
phase c =
  let (x, y) = unlift c
  in atan2 y x


-- | Return the real part of a complex number
--
real :: Elt a => Exp (Complex a) -> Exp a
real = A.fst

-- | Return the imaginary part of a complex number
--
imag :: Elt a => Exp (Complex a) -> Exp a
imag = A.snd

-- | Return the complex conjugate of a complex number, defined as
--
-- > conj(Z) = X - iY
--
conj :: (Elt a, IsNum a) => Exp (Complex a) -> Exp (Complex a)
conj z = lift (real z, - imag z)

