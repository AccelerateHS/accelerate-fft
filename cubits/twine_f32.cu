/*
 * Module      : Twine
 * Copyright   : [2016] Trevor L. McDonell
 * License     : BSD3
 *
 * Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
 * Stability   : experimental
 * Portability : non-portable (GHC extensions)
 *
 * Convert between Accelerate's Struct-of-Array representation of complex
 * numbers and the Array-of-Struct representation necessary for CUBLAS.
 *
 */

#include <cuda.h>
#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

__global__ void interleave
(
    cuFloatComplex * __restrict__ cplx,
    const float * __restrict__ real,
    const float * __restrict__ imag,
    const int size
)
{
    const int gridSize = blockDim.x * gridDim.x;
    int ix;

    for (ix = blockDim.x * blockIdx.x + threadIdx.x; ix < size; ix += gridSize) {
      const float re = real[ix];
      const float im = imag[ix];

      cplx[ix] = make_cuFloatComplex(re, im);
    }
}

__global__ void deinterleave
(
    float * __restrict__ real,
    float * __restrict__ imag,
    const cuFloatComplex * __restrict__ cplx,
    const int size
)
{
    const int gridSize = blockDim.x * gridDim.x;
    int ix;

    for (ix = blockDim.x * blockIdx.x + threadIdx.x; ix < size; ix += gridSize) {
      const cuFloatComplex c = cplx[ix];

      real[ix] = cuCrealf(c);
      imag[ix] = cuCimagf(c);
    }
}

#ifdef __cplusplus
}
#endif

