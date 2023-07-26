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

dd_real Rlansy(const char *norm, const char *uplo, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *work) {
    dd_real return_value = 0.0;
    //
    const dd_real zero = 0.0;
    dd_real value = 0.0;
    mplapackint j = 0;
    mplapackint i = 0;
    dd_real sum = 0.0;
    dd_real absa = 0.0;
    dd_real ssq[2];
    const dd_real one = 1.0;
    dd_real colssq[2];
    if (n == 0) {
        value = zero;
    } else if (Mlsame_dd(norm, "M")) {
        //
        //        Find std::max(abs(A(i,j))).
        //
        value = zero;
        if (Mlsame_dd(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    sum = abs(a[(i - 1) + (j - 1) * lda]);
                    if (value < sum || Risnan(sum)) {
                        value = sum;
                    }
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = j; i <= n; i = i + 1) {
                    sum = abs(a[(i - 1) + (j - 1) * lda]);
                    if (value < sum || Risnan(sum)) {
                        value = sum;
                    }
                }
            }
        }
    } else if ((Mlsame_dd(norm, "I")) || (Mlsame_dd(norm, "O")) || (Mlsame_dd(norm, "1"))) {
        //
        //        Find normI(A) ( = norm1(A), since A is symmetric).
        //
        value = zero;
        if (Mlsame_dd(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                sum = zero;
                for (i = 1; i <= j - 1; i = i + 1) {
                    absa = abs(a[(i - 1) + (j - 1) * lda]);
                    sum += absa;
                    work[i - 1] += absa;
                }
                work[j - 1] = sum + abs(a[(j - 1) + (j - 1) * lda]);
            }
            for (i = 1; i <= n; i = i + 1) {
                sum = work[i - 1];
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        } else {
            for (i = 1; i <= n; i = i + 1) {
                work[i - 1] = zero;
            }
            for (j = 1; j <= n; j = j + 1) {
                sum = work[j - 1] + abs(a[(j - 1) + (j - 1) * lda]);
                for (i = j + 1; i <= n; i = i + 1) {
                    absa = abs(a[(i - 1) + (j - 1) * lda]);
                    sum += absa;
                    work[i - 1] += absa;
                }
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        }
    } else if ((Mlsame_dd(norm, "F")) || (Mlsame_dd(norm, "E"))) {
        //
        //        Find normF(A).
        //        SSQ(1) is scale
        //        SSQ(2) is sum-of-squares
        //        For better accuracy, sum each column separately.
        //
        ssq[1 - 1] = zero;
        ssq[2 - 1] = one;
        //
        //        Sum off-diagonals
        //
        if (Mlsame_dd(uplo, "U")) {
            for (j = 2; j <= n; j = j + 1) {
                colssq[1 - 1] = zero;
                colssq[2 - 1] = one;
                Rlassq(j - 1, &a[(j - 1) * lda], 1, colssq[1 - 1], colssq[2 - 1]);
                Rcombssq(ssq, colssq);
            }
        } else {
            for (j = 1; j <= n - 1; j = j + 1) {
                colssq[1 - 1] = zero;
                colssq[2 - 1] = one;
                Rlassq(n - j, &a[((j + 1) - 1) + (j - 1) * lda], 1, colssq[1 - 1], colssq[2 - 1]);
                Rcombssq(ssq, colssq);
            }
        }
        ssq[2 - 1] = 2 * ssq[2 - 1];
        //
        //        Sum diagonal
        //
        colssq[1 - 1] = zero;
        colssq[2 - 1] = one;
        Rlassq(n, a, lda + 1, colssq[1 - 1], colssq[2 - 1]);
        Rcombssq(ssq, colssq);
        value = ssq[1 - 1] * sqrt(ssq[2 - 1]);
    }
    //
    return_value = value;
    return return_value;
    //
    //     End of Rlansy
    //
}
