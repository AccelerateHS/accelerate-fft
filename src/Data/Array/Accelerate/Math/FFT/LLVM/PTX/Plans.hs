{-# LANGUAGE MagicHash       #-}
{-# LANGUAGE RecordWildCards #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX.Plans
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX.Plans (

  Plans,
  createPlan,
  withPlan,

) where

import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.LLVM.PTX
import Data.Array.Accelerate.LLVM.PTX.Foreign

import Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base

import Control.Concurrent.MVar
import Control.Monad.State
import Data.HashMap.Strict
import qualified Data.HashMap.Strict                                as Map

import qualified Foreign.CUDA.Driver.Context                        as CUDA
import qualified Foreign.CUDA.FFT                                   as FFT

import GHC.Ptr
import GHC.Base
import Prelude                                                      hiding ( lookup )


data Plans a = Plans
  { plans   :: {-# UNPACK #-} !(MVar ( HashMap (Int, Int) (Lifetime FFT.Handle)))
  , create  :: a -> IO FFT.Handle
  , hash    :: a -> Int
  }


-- Create a new plan cache
--
{-# INLINE createPlan #-}
createPlan :: (a -> IO FFT.Handle) -> (a -> Int) -> IO (Plans a)
createPlan via mix =
  Plans <$> newMVar Map.empty <*> pure via <*> pure mix


-- Execute an operation with a cuFFT handle appropriate for the current
-- execution context.
--
-- Initial creation of the context is an atomic operation, but subsequently
-- multiple threads may use the context concurrently.
--
-- TLM: check that plans can be used concurrently
--
-- <http://docs.nvidia.com/cuda/cufft/index.html#thread-safety>
--
{-# INLINE withPlan #-}
withPlan :: Plans a -> a -> (FFT.Handle -> LLVM PTX b) -> LLVM PTX b
withPlan Plans{..} a k = do
  lc <- gets (deviceContext . ptxContext)
  h  <- liftIO $
          withLifetime lc  $ \ctx ->
          modifyMVar plans $ \pm  ->
            let key = (toKey ctx, hash a) in
            case Map.lookup key pm of
              -- handle does not exist yet; create it and add to the global
              -- state for reuse
              Nothing -> do
                h <- create a
                l <- newLifetime h
                addFinalizer lc $ modifyMVar plans (\pm' -> return (Map.delete key pm', ()))
                addFinalizer l  $ FFT.destroy h
                return ( Map.insert key l pm, l )

              -- return existing handle
              Just h  -> return (pm, h)
  --
  withLifetime' h k

{-# INLINE toKey #-}
toKey :: CUDA.Context -> Int
toKey (CUDA.Context (Ptr addr#)) = I# (addr2Int# addr#)

