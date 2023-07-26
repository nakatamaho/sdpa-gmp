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

#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

void Rorg2r(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint &info) {
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (k < 0 || k > n) {
        info = -3;
    } else if (lda < std::max((mplapackint)1, m)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla_gmp("Rorg2r", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Initialise columns k+1:n to columns of the unit matrix
    //
    mplapackint j = 0;
    mplapackint l = 0;
    const mpf_class zero = 0.0;
    const mpf_class one = 1.0;
    for (j = k + 1; j <= n; j = j + 1) {
        for (l = 1; l <= m; l = l + 1) {
            a[(l - 1) + (j - 1) * lda] = zero;
        }
        a[(j - 1) + (j - 1) * lda] = one;
    }
    //
    mplapackint i = 0;
    for (i = k; i >= 1; i = i - 1) {
        //
        //        Apply H(i) to A(i:m,i:n) from the left
        //
        if (i < n) {
            a[(i - 1) + (i - 1) * lda] = one;
            Rlarf("Left", m - i + 1, n - i, &a[(i - 1) + (i - 1) * lda], 1, tau[i - 1], &a[(i - 1) + ((i + 1) - 1) * lda], lda, work);
        }
        if (i < m) {
            Rscal(m - i, -tau[i - 1], &a[((i + 1) - 1) + (i - 1) * lda], 1);
        }
        a[(i - 1) + (i - 1) * lda] = one - tau[i - 1];
        //
        //        Set A(1:i-1,i) to zero
        //
        for (l = 1; l <= i - 1; l = l + 1) {
            a[(l - 1) + (i - 1) * lda] = zero;
        }
    }
    //
    //     End of Rorg2r
    //
}
