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

void Rsytd2(const char *uplo, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *d, dd_real *e, dd_real *tau, mplapackint &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    bool upper = Mlsame_dd(uplo, "U");
    if (!upper && !Mlsame_dd(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla_dd("Rsytd2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    mplapackint i = 0;
    dd_real taui = 0.0;
    const dd_real zero = 0.0;
    const dd_real one = 1.0;
    const dd_real half = 1.0 / 2.0;
    dd_real alpha = 0.0;
    if (upper) {
        //
        //        Reduce the upper triangle of A
        //
        for (i = n - 1; i >= 1; i = i - 1) {
            //
            //           Generate elementary reflector H(i) = I - tau * v * v**T
            //           to annihilate A(1:i-1,i+1)
            //
            Rlarfg(i, a[(i - 1) + ((i + 1) - 1) * lda], &a[((i + 1) - 1) * lda], 1, taui);
            e[i - 1] = a[(i - 1) + ((i + 1) - 1) * lda];
            //
            if (taui != zero) {
                //
                //              Apply H(i) from both sides to A(1:i,1:i)
                //
                a[(i - 1) + ((i + 1) - 1) * lda] = one;
                //
                //              Compute  x := tau * A * v  storing x in TAU(1:i)
                //
                Rsymv(uplo, i, taui, a, lda, &a[((i + 1) - 1) * lda], 1, zero, tau, 1);
                //
                //              Compute  w := x - 1/2 * tau * (x**T * v) * v
                //
                alpha = -half * taui * Rdot(i, tau, 1, &a[((i + 1) - 1) * lda], 1);
                Raxpy(i, alpha, &a[((i + 1) - 1) * lda], 1, tau, 1);
                //
                //              Apply the transformation as a rank-2 update:
                //                 A := A - v * w**T - w * v**T
                //
                Rsyr2(uplo, i, -one, &a[((i + 1) - 1) * lda], 1, tau, 1, a, lda);
                //
                a[(i - 1) + ((i + 1) - 1) * lda] = e[i - 1];
            }
            d[(i + 1) - 1] = a[((i + 1) - 1) + ((i + 1) - 1) * lda];
            tau[i - 1] = taui;
        }
        d[1 - 1] = a[(1 - 1)];
    } else {
        //
        //        Reduce the lower triangle of A
        //
        for (i = 1; i <= n - 1; i = i + 1) {
            //
            //           Generate elementary reflector H(i) = I - tau * v * v**T
            //           to annihilate A(i+2:n,i)
            //
	    Rlarfg(n - i, a[((i + 1) - 1) + (i - 1) * lda], &a[(std::min(i + 2, n) - 1) + (i - 1) * lda], 1, taui);
            e[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
            //
            if (taui != zero) {
                //
                //              Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                a[((i + 1) - 1) + (i - 1) * lda] = one;
                //
                //              Compute  x := tau * A * v  storing y in TAU(i:n-1)
                //
                Rsymv(uplo, n - i, taui, &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda, &a[((i + 1) - 1) + (i - 1) * lda], 1, zero, &tau[i - 1], 1);
                //
                //              Compute  w := x - 1/2 * tau * (x**T * v) * v
                //
                alpha = -half * taui * Rdot(n - i, &tau[i - 1], 1, &a[((i + 1) - 1) + (i - 1) * lda], 1);
                Raxpy(n - i, alpha, &a[((i + 1) - 1) + (i - 1) * lda], 1, &tau[i - 1], 1);
                //
                //              Apply the transformation as a rank-2 update:
                //                 A := A - v * w**T - w * v**T
                //
                Rsyr2(uplo, n - i, -one, &a[((i + 1) - 1) + (i - 1) * lda], 1, &tau[i - 1], 1, &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda);
                //
                a[((i + 1) - 1) + (i - 1) * lda] = e[i - 1];
            }
            d[i - 1] = a[(i - 1) + (i - 1) * lda];
            tau[i - 1] = taui;
        }
        d[n - 1] = a[(n - 1) + (n - 1) * lda];
    }
    //
    //     End of Rsytd2
    //
}
