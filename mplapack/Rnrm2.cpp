/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#include <mpblas_dd.h>

dd_real Rnrm2(mplapackint const n, dd_real *x, mplapackint const incx) {
    dd_real return_value = 0.0;
    //
    //  -- Reference BLAS level1 routine --
    //  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    const dd_real zero = 0.0;
    dd_real norm = 0.0;
    dd_real scale = 0.0;
    const dd_real one = 1.0;
    dd_real ssq = 0.0;
    mplapackint ix = 0;
    dd_real absxi = 0.0;
    if (n < 1 || incx < 1) {
        norm = zero;
    } else if (n == 1) {
        norm = abs(x[1 - 1]);
    } else {
        scale = zero;
        ssq = one;
        //        The following loop is equivalent to this call to the LAPACK
        //        auxiliary routine:
        //        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
        //
        for (ix = 1; ix <= 1 + (n - 1) * incx; ix = ix + incx) {
            if (x[ix - 1] != zero) {
                absxi = abs(x[ix - 1]);
                if (scale < absxi) {
                    ssq = one + ssq * pow2((scale / absxi));
                    scale = absxi;
                } else {
                    ssq += pow2((absxi / scale));
                }
            }
        }
        norm = scale * sqrt(ssq);
    }
    //
    return_value = norm;
    return return_value;
    //
    //     End of Rnrm2.
    //
}
