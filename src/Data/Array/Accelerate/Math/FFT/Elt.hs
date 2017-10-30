{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.Elt
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.Elt
  where

import Data.Array.Accelerate                                        as A
import Prelude                                                      as P


-- | The type of supported FFT elements; namely 'Float' and 'Double'.
--
type FFTElt e = (P.Num e, A.RealFloat e, A.FromIntegral Int e, A.IsFloating e)

