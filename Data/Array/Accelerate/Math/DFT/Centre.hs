{-# LANGUAGE TypeOperators #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT.Centre
-- Copyright   : [2012..2013] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell, Robert Clifton-Everest
-- License     : BSD3
--
-- Maintainer  : Manuel M T Chakravarty <chak@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- These transforms allow the centering of the frequency domain of a DFT such
-- that the the zero frequency is in the middle. The centering transform, when
-- performed on the input of a DFT, will cause zero frequency to be centred in
-- the middle. The shifting transform however takes the output of a DFT to
-- give the same result. Therefore the relationship between the two is:
--
-- > fft(center(X)) = shift(fft(X))
--
module Data.Array.Accelerate.Math.DFT.Centre (

  centre1D, centre2D, centre3D,
  shift1D,  shift2D,  shift3D,

) where

import Prelude                                  as P
import Data.Array.Accelerate                    as A
import Data.Array.Accelerate.Math.Complex


-- | Apply the centring transform to a vector
--
centre1D :: (Elt e, IsFloating e)
         => Acc (Array DIM1 (Complex e))
         -> Acc (Array DIM1 (Complex e))
centre1D arr
  = A.generate (shape arr)
               (\ix -> let Z :. x = unlift ix           :: Z :. Exp Int
                       in  lift (-1 ** A.fromIntegral x, A.constant 0) * arr!ix)

-- | Apply the centring transform to a matrix
--
centre2D :: (Elt e, IsFloating e)
         => Acc (Array DIM2 (Complex e))
         -> Acc (Array DIM2 (Complex e))
centre2D arr
  = A.generate (shape arr)
               (\ix -> let Z :. y :. x = unlift ix      :: Z :. Exp Int :. Exp Int
                       in  lift (-1 ** A.fromIntegral (y + x), A.constant 0) * arr!ix)

-- | Apply the centring transform to a 3D array
--
centre3D :: (Elt e, IsFloating e)
         => Acc (Array DIM3 (Complex e))
         -> Acc (Array DIM3 (Complex e))
centre3D arr
  = A.generate (shape arr)
               (\ix -> let Z :. z :. y :. x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
                       in  lift (-1 ** A.fromIntegral (z + y + x), A.constant 0) * arr!ix)


-- | Apply the shifting transform to a vector
--
shift1D :: Elt e => Acc (Vector e) -> Acc (Vector e)
shift1D arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = let Z:.x = unlift ix :: Z :. Exp Int
        in index1 (x <* mw ? (x + mw, x - mw))
    Z:.w    = unlift (A.shape arr)
    mw      = w `div` 2


-- | Apply the shifting transform to a 2D array
--
shift2D :: Elt e => Acc (Array DIM2 e) -> Acc (Array DIM2 e)
shift2D arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = let Z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int
        in index2 (y <* mh ? (y + mh, y - mh))
                  (x <* mw ? (x + mw, x - mw))
    Z:.h:.w = unlift (A.shape arr)
    (mh,mw) = (h `div` 2, w `div` 2)


-- | Apply the shifting transform to a 3D array
--
shift3D :: Elt e => Acc (Array DIM3 e) -> Acc (Array DIM3 e)
shift3D arr
  = A.backpermute (A.shape arr) p arr
  where
    p ix
      = let Z:.z:.y:.x = unlift ix :: Z :. Exp Int :. Exp Int :. Exp Int
        in index3 (z <* md ? (z + md, z - md))
                  (y <* mh ? (y + mh, y - mh))
                  (x <* mw ? (x + mw, x - mw))
    Z:.h:.w:.d = unlift (A.shape arr)
    (mh,mw,md) = (h `div` 2, w `div` 2, d `div` 2)
    index3 i j k = lift (Z:.i:.j:.k)

