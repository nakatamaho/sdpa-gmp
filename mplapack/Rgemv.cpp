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

void Rgemv(const char *trans, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy) {
    //
    //  -- Reference BLAS level2 routine --
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //
    //     Test the input parameters.
    //
    mplapackint info = 0;
    if (!Mlsame_dd(trans, "N") && !Mlsame_dd(trans, "T") && !Mlsame_dd(trans, "C")) {
        info = 1;
    } else if (m < 0) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (lda < std::max((mplapackint)1, m)) {
        info = 6;
    } else if (incx == 0) {
        info = 8;
    } else if (incy == 0) {
        info = 11;
    }
    if (info != 0) {
        Mxerbla_dd("Rgemv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const dd_real zero = 0.0;
    const dd_real one = 1.0;
    if ((m == 0) || (n == 0) || ((alpha == zero) && (beta == one))) {
        return;
    }
    //
    //     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    //     up the start points in  X  and  Y.
    //
    mplapackint lenx = 0;
    mplapackint leny = 0;
    if (Mlsame_dd(trans, "N")) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    mplapackint kx = 0;
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * incx;
    }
    mplapackint ky = 0;
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * incy;
    }
    //
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through A.
    //
    //     First form  y := beta*y.
    //
    mplapackint i = 0;
    mplapackint iy = 0;
    if (beta != one) {
        if (incy == 1) {
            if (beta == zero) {
                for (i = 1; i <= leny; i = i + 1) {
                    y[i - 1] = zero;
                }
            } else {
                for (i = 1; i <= leny; i = i + 1) {
                    y[i - 1] = beta * y[i - 1];
                }
            }
        } else {
            iy = ky;
            if (beta == zero) {
                for (i = 1; i <= leny; i = i + 1) {
                    y[iy - 1] = zero;
                    iy += incy;
                }
            } else {
                for (i = 1; i <= leny; i = i + 1) {
                    y[iy - 1] = beta * y[iy - 1];
                    iy += incy;
                }
            }
        }
    }
    if (alpha == zero) {
        return;
    }
    mplapackint jx = 0;
    mplapackint j = 0;
    dd_real temp = 0.0;
    mplapackint jy = 0;
    mplapackint ix = 0;
    if (Mlsame_dd(trans, "N")) {
        //
        //        Form  y := alpha*A*x + y.
        //
        jx = kx;
        if (incy == 1) {
            for (j = 1; j <= n; j = j + 1) {
                temp = alpha * x[jx - 1];
                for (i = 1; i <= m; i = i + 1) {
                    y[i - 1] += temp * a[(i - 1) + (j - 1) * lda];
                }
                jx += incx;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                temp = alpha * x[jx - 1];
                iy = ky;
                for (i = 1; i <= m; i = i + 1) {
                    y[iy - 1] += temp * a[(i - 1) + (j - 1) * lda];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        //
        //        Form  y := alpha*A**T*x + y.
        //
        jy = ky;
        if (incx == 1) {
            for (j = 1; j <= n; j = j + 1) {
                temp = zero;
                for (i = 1; i <= m; i = i + 1) {
                    temp += a[(i - 1) + (j - 1) * lda] * x[i - 1];
                }
                y[jy - 1] += alpha * temp;
                jy += incy;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                temp = zero;
                ix = kx;
                for (i = 1; i <= m; i = i + 1) {
                    temp += a[(i - 1) + (j - 1) * lda] * x[ix - 1];
                    ix += incx;
                }
                y[jy - 1] += alpha * temp;
                jy += incy;
            }
        }
    }
    //
    //     End of Rgemv .
    //
}
