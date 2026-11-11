{-# LANGUAGE MagicHash       #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TupleSections #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX.Plans
-- Copyright   : [2017..2020] The Accelerate Team
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX.Plans (

  Plans,
  createPlan,
  withPlan,

) where

import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.LLVM.PTX hiding (stream, poll)
import Data.Array.Accelerate.LLVM.PTX.Foreign

import Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base

import Control.Concurrent.MVar
import Control.Monad.Catch
import Control.Monad.State
import Data.HashMap.Strict hiding (map, update)
import qualified Data.HashMap.Strict                                as Map

import qualified Foreign.CUDA.Driver.Context                        as CUDA
import qualified Foreign.CUDA.Driver.Stream                         as CUDA
import qualified Foreign.CUDA.FFT                                   as FFT

import GHC.Ptr
import GHC.Base
import Prelude                                                      hiding ( lookup, mapM )
import Data.Maybe
import Control.Arrow (second)
import Data.Function ((&))
import Control.Monad.Reader (asks)


data Plans a = Plans
  { plans   :: {-# UNPACK #-} !(MVar ( HashMap (Int, Int) [(Lifetime FFT.Handle, Maybe (Par PTX Bool, CUDA.Stream))]))
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
-- TODO: Determine if this handle is used in the same stream.
{-# INLINE withPlan #-}
withPlan :: Plans a -> a -> (FFT.Handle -> Par PTX (Future b)) -> Par PTX (Future b)
withPlan Plans{..} a k = do
  lc <- gets (deviceContext . ptxContext)
  ls <- asks ptxStream
  withLifetime' ls $ \stream ->
    withLifetime' lc  $ \ctx -> do
      let key = (toKey ctx, hash a)
      -- Extract an existing cuFFT plan handle from our plan cache that isn't busy,
      -- if one cannot be found, create a new cuFFT handle.
      h <- modifyMVar' plans $ \pm -> do
              let maybeHandles = pm !? key
                  handles = fromMaybe [] maybeHandles

                  update Nothing = pure Nothing
                  update orig@(Just (isReady, _)) = isReady >>= \case
                    True  -> pure Nothing
                    False -> pure orig

              updatedHandles <- zip (map fst handles) <$> mapM (update . snd) handles

              -- Extract first handle which is either entirely ready or is used but within the same stream
              let extractFirstReady []                                    = (Nothing, [])
                  extractFirstReady (x@(_, Nothing):xs)                   = (Just x, xs)
                  extractFirstReady (x@(_, Just (_, s)):xs) | stream == s = (Just x, xs)
                  extractFirstReady (x@(_, Just _):xs)                    = second (x:) $ extractFirstReady xs

                  (maybeReadyHandle, otherHandles) = extractFirstReady updatedHandles

                  newHandle = liftIO $ do
                    h <- create a
                    l <- newLifetime h
                    addFinalizer l $ FFT.destroy h
                    when (isNothing maybeHandles) $
                      addFinalizer lc $ modifyMVar_ plans $ pure . Map.delete key
                    pure l

              maybeReadyHandle & maybe newHandle (pure . fst)
                               & fmap (Map.insert key otherHandles pm,)
      -- Ensure the handle is always returned back to the plan cache
      let returnHandle = liftIO $ modifyMVar_ plans $ pure . Map.adjust ((h, Nothing):) key
      flip onException returnHandle $ do
        -- Invoke user-provided function with cuFFT handle
        future <- withLifetime' h k
        -- Push new cuFFT plan-handle onto list of plan-handles of equal settings,
        -- w/ callback to check if the cuFFT handle is ready to use again.
        planHandleEntry <- (h,) . Just . (,stream) . fmap isJust . poll <$> statusHandle future
        liftIO $ modifyMVar_ plans $ pure . Map.adjust (planHandleEntry:) key
        pure future

{-# INLINE toKey #-}
toKey :: CUDA.Context -> Int
toKey (CUDA.Context (Ptr addr#)) = I# (addr2Int# addr#)
