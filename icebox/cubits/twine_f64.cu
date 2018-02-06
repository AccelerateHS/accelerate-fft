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
    cuDoubleComplex * __restrict__ cplx,
    const double * __restrict__ real,
    const double * __restrict__ imag,
    const int size
)
{
    const int gridSize = blockDim.x * gridDim.x;
    int ix;

    for (ix = blockDim.x * blockIdx.x + threadIdx.x; ix < size; ix += gridSize) {
      const double re = real[ix];
      const double im = imag[ix];

      cplx[ix] = make_cuDoubleComplex(re, im);
    }
}

__global__ void deinterleave
(
    double * __restrict__ real,
    double * __restrict__ imag,
    const cuDoubleComplex * __restrict__ cplx,
    const int size
)
{
    const int gridSize = blockDim.x * gridDim.x;
    int ix;

    for (ix = blockDim.x * blockIdx.x + threadIdx.x; ix < size; ix += gridSize) {
      const cuDoubleComplex c = cplx[ix];

      real[ix] = cuCreal(c);
      imag[ix] = cuCimag(c);
    }
}

#ifdef __cplusplus
}
#endif

