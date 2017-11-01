{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE PolyKinds           #-}
{-# LANGUAGE ScopedTypeVariables #-}
-- |
-- Module      : Test.ShowType
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Test.ShowType
  where

import Data.Complex

data ArgType (a :: *) = AT

showType :: forall proxy a. Show (ArgType a) => proxy a -> String
showType _ = show (AT :: ArgType a)

instance Show (ArgType a) => Show (ArgType (Complex a)) where
  show _ = "Complex " ++ show (AT :: ArgType a)

instance Show (ArgType Float)  where show _ = "Float"
instance Show (ArgType Double) where show _ = "Double"

