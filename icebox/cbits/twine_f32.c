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
 * numbers and the Array-of-Struct representation used by BLAS.
 */

#include <complex.h>
#include "HsFFI.h"

#ifdef __cplusplus
extern "C" {
#endif

void interleave_f32
(
    const StgInt start,
    const StgInt end,
    complex float * __restrict__ cplx,
    const float * __restrict__ real,
    const float * __restrict__ imag
)
{
    StgInt i;
    for (i = start; i < end; ++i) {
        const float re = real[i];
        const float im = imag[i];

        cplx[i] = re + im * I;
    }
}

void deinterleave_f32
(
    const StgInt start,
    const StgInt end,
    float * __restrict__ real,
    float * __restrict__ imag,
    const complex float * __restrict__ cplx
)
{
    StgInt i;
    for (i = start; i < end; ++i) {
        const complex float c = cplx[i];

        real[i] = crealf(c);
        imag[i] = cimagf(c);
    }
}

#ifdef __cplusplus
}
#endif

