{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE MagicHash           #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
-- |
-- Module      : Data.Array.Accelerate.Math.FFT.LLVM.PTX.Twine
-- Copyright   : [2017] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <trevor.mcdonell@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Data.Array.Accelerate.Math.FFT.LLVM.PTX.Twine (

  interleave,
  deinterleave,

) where

import Data.Array.Accelerate.Array.Sugar
import Data.Array.Accelerate.Data.Complex
import Data.Array.Accelerate.Lifetime
import Data.Array.Accelerate.Type

import Data.Array.Accelerate.LLVM.PTX
import Data.Array.Accelerate.LLVM.PTX.Foreign

import Data.Array.Accelerate.Math.FFT.LLVM.PTX.Base

import Control.Concurrent.MVar
import Control.Monad.State
import Data.ByteString                                              ( ByteString )
import Data.FileEmbed
import Data.IntMap.Strict                                           ( IntMap )
import System.IO.Unsafe
import qualified Data.IntMap.Strict                                 as IM

import Foreign.CUDA.Ptr                                             ( DevicePtr )
import Foreign.CUDA.Analysis
import qualified Foreign.CUDA.Driver                                as CUDA
import qualified Foreign.CUDA.Driver.Stream                         as CUDA

import GHC.Ptr
import GHC.Base
import Prelude                                                      hiding ( lookup )


interleave
    :: forall sh e a b. (Shape sh, Elt e, IsFloating e, DevicePtrs e ~ DevicePtr a)
    => Array sh (Complex e)
    -> Stream
    -> (DevicePtr a -> LLVM PTX b)  -- device pointer is in packed representation
    -> LLVM PTX b
interleave arr s k = do
  let n = size (shape arr)
  --
  cplx <- allocateRemote (Z :. n*2) :: LLVM PTX (Vector e)  -- packed representation
  withTwine (floatingType :: FloatingType e)  $ \(_,pack,_) -> do
    withComplexArrayPtrs arr s                $ \d_re d_im  -> do
      withScalarArrayPtr cplx s               $ \d_cplx     -> do
        withLifetime' s                       $ \s'         -> do
          liftIO $ launch pack s' n d_cplx d_re d_im
        k (CUDA.castDevPtr d_cplx)


deinterleave
    :: forall sh e a. (Shape sh, IsFloating e, DevicePtrs e ~ DevicePtr a)
    => Array sh (Complex e)   -- [out] destination
    -> DevicePtr a            -- [in]  in packed representation
    -> Stream
    -> LLVM PTX ()
deinterleave arr d_cplx s = do
  let n = size (shape arr)
  --
  withTwine (floatingType :: FloatingType e)  $ \(_,_,unpack) -> do
    withComplexArrayPtrs arr s                $ \d_re d_im    -> do
      withLifetime' s                         $ \s'           -> do
        liftIO $ launch unpack s' n d_re d_im d_cplx


withTwine
    :: FloatingType e
    -> ((CUDA.Module, Kernel, Kernel) -> LLVM PTX b)
    -> LLVM PTX b
withTwine tR k = do
  ptx <- gets ptxContext
  let lc  = deviceContext ptx
      prp = deviceProperties ptx
      mds = modules tR
  --
  mdl <- liftIO $ do
    withLifetime lc $ \ctx -> do
     modifyMVar mds $ \im  -> do
      let key = toKey ctx
      case IM.lookup key im of
        -- Module is not loaded yet; add to the current context and the global
        -- state for later reuse
        Nothing -> do
          mdl     <- CUDA.loadData $ case tR of
                                       TypeFloat{}   -> ptx_twine_f32
                                       TypeDouble{}  -> ptx_twine_f64
                                       TypeCFloat{}  -> ptx_twine_f32
                                       TypeCDouble{} -> ptx_twine_f64
          pack    <- mkKernel "interleave"   mdl prp
          unpack  <- mkKernel "deinterleave" mdl prp
          let mkk = (mdl, pack, unpack)
          --
          lm      <- newLifetime mkk
          addFinalizer lc $ modifyMVar mds (\im' -> return (IM.delete key im', ()))
          addFinalizer lm $ CUDA.unload mdl
          return ( IM.insert key lm im, lm )

        -- Return existing module
        Just lm  -> return (im, lm)
  --
  withLifetime' mdl k


toKey :: CUDA.Context -> IM.Key
toKey (CUDA.Context (Ptr addr#)) = I# (addr2Int# addr#)


launch :: Kernel -> CUDA.Stream -> Int -> DevicePtr e -> DevicePtr e -> DevicePtr e -> IO ()
launch Kernel{..} s n dx dy dz =
  CUDA.launchKernel kernelFun (kernelThreadBlocks n,1,1) (kernelThreadBlockSize,1,1) kernelSharedMemBytes (Just s)
    [ CUDA.VArg dx, CUDA.VArg dy, CUDA.VArg dz, CUDA.IArg (fromIntegral n) ]

mkKernel :: String -> CUDA.Module -> CUDA.DeviceProperties -> IO Kernel
mkKernel name mdl prp = do
  fun <- CUDA.getFun mdl name
  reg <- CUDA.requires fun CUDA.NumRegs
  let
      blockSize   = 256
      sharedMem   = 0
      maxBlocks   = maxResidentBlocks prp blockSize reg sharedMem
      numBlocks n = maxBlocks `min` ((n + blockSize - 1) `quot` blockSize)
  --
  return $ Kernel fun sharedMem blockSize numBlocks name

data Kernel = Kernel {
    kernelFun               :: {-# UNPACK #-} !CUDA.Fun
  , kernelSharedMemBytes    :: {-# UNPACK #-} !Int
  , kernelThreadBlockSize   :: {-# UNPACK #-} !Int
  , kernelThreadBlocks      :: (Int -> Int)
  , kernelName              :: String
  }

modules :: FloatingType e -> MVar (IntMap (Lifetime (CUDA.Module, Kernel, Kernel)))
modules TypeFloat{}   = modules_f32
modules TypeDouble{}  = modules_f64
modules TypeCFloat{}  = modules_f32
modules TypeCDouble{} = modules_f64

{-# NOINLINE modules_f32 #-}
modules_f32 :: MVar (IntMap (Lifetime (CUDA.Module, Kernel, Kernel)))
modules_f32 = unsafePerformIO $ newMVar IM.empty

{-# NOINLINE modules_f64 #-}
modules_f64 :: MVar (IntMap (Lifetime (CUDA.Module, Kernel, Kernel)))
modules_f64 = unsafePerformIO $ newMVar IM.empty

ptx_twine_f32 :: ByteString
ptx_twine_f32 = $(makeRelativeToProject "cubits/twine_f32.ptx" >>= embedFile)

ptx_twine_f64 :: ByteString
ptx_twine_f64 = $(makeRelativeToProject "cubits/twine_f64.ptx" >>= embedFile)

