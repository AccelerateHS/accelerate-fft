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

void interleave_f64
(
    const StgInt start,
    const StgInt end,
    complex double * __restrict__ cplx,
    const double * __restrict__ real,
    const double * __restrict__ imag
)
{
    StgInt i;
    for (i = start; i < end; ++i) {
        const double re = real[i];
        const double im = imag[i];

        cplx[i] = re + im * I;
    }
}

void deinterleave_f64
(
    const StgInt start,
    const StgInt end,
    double * __restrict__ real,
    double * __restrict__ imag,
    const complex double * __restrict__ cplx
)
{
    StgInt i;
    for (i = start; i < end; ++i) {
        const complex double c = cplx[i];

        real[i] = creal(c);
        imag[i] = cimag(c);
    }
}

#ifdef __cplusplus
}
#endif

