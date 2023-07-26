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
#include <mplapack_dd.h>

void Rpotf2(const char *uplo, mplapackint const n, dd_real *a, mplapackint const lda, mplapackint &info) {
    bool upper = false;
    mplapackint j = 0;
    dd_real ajj = 0.0;
    const dd_real zero = 0.0;
    const dd_real one = 1.0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    upper = Mlsame_dd(uplo, "U");
    if (!upper && !Mlsame_dd(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla_dd("Rpotf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    if (upper) {
        //
        //        Compute the Cholesky factorization A = U**T *U.
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //           Compute U(J,J) and test for non-positive-definiteness.
            //
            ajj = a[(j - 1) + (j - 1) * lda] - Rdot(j - 1, &a[(j - 1) * lda], 1, &a[(j - 1) * lda], 1);
            if (ajj <= zero || Risnan(ajj)) {
                a[(j - 1) + (j - 1) * lda] = ajj;
                goto statement_30;
            }
            ajj = sqrt(ajj);
            a[(j - 1) + (j - 1) * lda] = ajj;
            //
            //           Compute elements J+1:N of row J.
            //
            if (j < n) {
                Rgemv("Transpose", j - 1, n - j, -one, &a[((j + 1) - 1) * lda], lda, &a[(j - 1) * lda], 1, one, &a[(j - 1) + ((j + 1) - 1) * lda], lda);
                Rscal(n - j, one / ajj, &a[(j - 1) + ((j + 1) - 1) * lda], lda);
            }
        }
    } else {
        //
        //        Compute the Cholesky factorization A = L*L**T.
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //           Compute L(J,J) and test for non-positive-definiteness.
            //
            ajj = a[(j - 1) + (j - 1) * lda] - Rdot(j - 1, &a[(j - 1)], lda, &a[(j - 1)], lda);
            if (ajj <= zero || Risnan(ajj)) {
                a[(j - 1) + (j - 1) * lda] = ajj;
                goto statement_30;
            }
            ajj = sqrt(ajj);
            a[(j - 1) + (j - 1) * lda] = ajj;
            //
            //           Compute elements J+1:N of column J.
            //
            if (j < n) {
                Rgemv("No transpose", n - j, j - 1, -one, &a[((j + 1) - 1)], lda, &a[(j - 1)], lda, one, &a[((j + 1) - 1) + (j - 1) * lda], 1);
                Rscal(n - j, one / ajj, &a[((j + 1) - 1) + (j - 1) * lda], 1);
            }
        }
    }
    goto statement_40;
//
statement_30:
    info = j;
//
statement_40:;
    //
    //     End of Rpotf2
    //
}
