{-# LANGUAGE TypeOperators #-}
-- |
-- Module      : Data.Array.Accelerate.Math.DFT.Centre
-- Copyright   : [2012] Manuel M T Chakravarty, Gabriele Keller, Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Manuel M T Chakravarty <chak@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Applying these transforms to the input of a DFT causes the output to be
-- centred so that the zero frequency is in the middle.
--
module Data.Array.Accelerate.Math.DFT.Centre (

  centre1D, centre2D, centre3D,

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

