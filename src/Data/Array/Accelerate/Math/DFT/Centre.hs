{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT.Centre
-- Copyright   : [2012..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- These transforms allow the centering of the frequency domain of a DFT such
-- that the zero frequency is in the middle. The centering transform, when
-- performed on the input of a DFT, will cause zero frequency to be centred in
-- the middle. The shifting transform however takes the output of a DFT to give
-- the same result. Therefore the relationship between the two is:
--
-- > fft(center(X)) = shift(fft(X))
--
module Data.Array.Accelerate.Math.DFT.Centre (

  centre1D, centre2D, centre3D,
  shift1D,  shift2D,  shift3D,
  ishift1D,  ishift2D,  ishift3D,

) where

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Data.Complex


-- | Apply the centring transform to a vector
--
centre1D
    :: (A.RealFloat e, A.FromIntegral Int e)
    => Acc (Array DIM1 (Complex e))
    -> Acc (Array DIM1 (Complex e))
centre1D arr
  = A.generate (shape arr)
               (\ix -> let Z :. x = unlift ix           :: Z :. Exp Int
                       in  lift (((-1) ** A.fromIntegral x) :+ 0) * arr!ix)

-- | Apply the centring transform to a matrix
--
centre2D
    :: (A.RealFloat e, A.FromIntegral Int e)
    => Acc (Array DIM2 (Complex e))
    -> Acc (Array DIM2 (Complex e))
centre2D arr
  = A.generate (shape arr)
               (\ix -> let Z :. y :. x = unlift ix      :: Z :. Exp Int :. Exp Int
                       in  lift (((-1) ** A.fromIntegral (y + x)) :+ 0) * arr!ix)

-- | Apply the centring transform to a 3D array
--
centre3D
    :: (A.RealFloat e, A.FromIntegral Int e)
    => Acc (Array DIM3 (Complex e))
    -> Acc (Array DIM3 (Complex e))
centre3D arr
  = A.generate (shape arr)
               (\ix -> let Z :. z :. y :. x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
                       in  lift (((-1) ** A.fromIntegral (z + y + x)) :+ 0) * arr!ix)


-- | Apply the shifting transform to a vector
--
shift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D arr = backpermute sh p arr
      where
        sh      = shape arr
        n       = indexHead sh
        --
        shift   = (n `quot` 2) + boolToInt (A.odd n)
        roll i  = (i+shift) `rem` n
        p       = ilift1 roll

-- | The inverse of the shift1D function, such that
-- > ishift1D (shift1D v) = ishift1D (shift1D v) = v
-- for all vectors
--
ishift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
ishift1D arr = backpermute sh p arr
      where
        sh      = shape arr
        n       = indexHead sh
        --
        shift   = (n `quot` 2)-- + boolToInt (A.odd n)
        roll i  = (i+shift) `rem` n
        p       = ilift1 roll

-- | Apply the shifting transform to a 2D array
--
shift2D :: Elt e => Acc (Array DIM2 e) -> Acc (Array DIM2 e)
shift2D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. h :. w = unlift sh
    --
    shifth = (h `quot` 2) + boolToInt (A.odd h)
    shiftw = (w `quot` 2) + boolToInt (A.odd w)

    p ix
      = let Z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int
        in index2 ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

-- | The inverse of the shift2D function
--
ishift2D :: Elt e => Acc (Array DIM2 e) -> Acc (Array DIM2 e)
ishift2D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. h :. w = unlift sh
    --
    shifth = (h `quot` 2)
    shiftw = (w `quot` 2)

    p ix
      = let Z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int
        in index2 ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

-- | Apply the shifting transform to a 3D array
--
shift3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
shift3D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. d :. h :. w = unlift sh
    --
    shiftd = (d `quot` 2) + boolToInt (A.odd d)
    shifth = (h `quot` 2) + boolToInt (A.odd h)
    shiftw = (w `quot` 2) + boolToInt (A.odd w)

    p ix
      = let Z:.z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
        in index3 ((z + shiftd) `rem` d)
                  ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)

-- | The inverse of the shift3D function
--
ishift3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
ishift3D arr
  = backpermute sh p arr
  where
    sh      = shape arr
    Z :. d :. h :. w = unlift sh
    --
    shiftd = (d `quot` 2)
    shifth = (h `quot` 2)
    shiftw = (w `quot` 2)

    p ix
      = let Z:.z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
        in index3 ((z + shiftd) `rem` d)
                  ((y + shifth) `rem` h)
                  ((x + shiftw) `rem` w)
