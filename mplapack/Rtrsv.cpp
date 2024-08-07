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

void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx) {
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
    if (!Mlsame_gmp(uplo, "U") && !Mlsame_gmp(uplo, "L")) {
        info = 1;
    } else if (!Mlsame_gmp(trans, "N") && !Mlsame_gmp(trans, "T") && !Mlsame_gmp(trans, "C")) {
        info = 2;
    } else if (!Mlsame_gmp(diag, "U") && !Mlsame_gmp(diag, "N")) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = 6;
    } else if (incx == 0) {
        info = 8;
    }
    if (info != 0) {
        Mxerbla_gmp("Rtrsv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    bool nounit = Mlsame_gmp(diag, "N");
    //
    //     Set up the start point in X if the increment is not unity. This
    //     will be  ( N - 1 )*INCX  too small for descending loops.
    //
    mplapackint kx = 0;
    if (incx <= 0) {
        kx = 1 - (n - 1) * incx;
    } else if (incx != 1) {
        kx = 1;
    }
    //
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through A.
    //
    mplapackint j = 0;
    const mpf_class zero = 0.0;
    mpf_class temp = 0.0;
    mplapackint i = 0;
    mplapackint jx = 0;
    mplapackint ix = 0;
    if (Mlsame_gmp(trans, "N")) {
        //
        //        Form  x := inv( A )*x.
        //
        if (Mlsame_gmp(uplo, "U")) {
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    if (x[j - 1] != zero) {
                        if (nounit) {
                            x[j - 1] = x[j - 1] / a[(j - 1) + (j - 1) * lda];
                        }
                        temp = x[j - 1];
                        for (i = j - 1; i >= 1; i = i - 1) {
                            x[i - 1] = x[i - 1] - temp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            } else {
                jx = kx + (n - 1) * incx;
                for (j = n; j >= 1; j = j - 1) {
                    if (x[jx - 1] != zero) {
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] / a[(j - 1) + (j - 1) * lda];
                        }
                        temp = x[jx - 1];
                        ix = jx;
                        for (i = j - 1; i >= 1; i = i - 1) {
                            ix = ix - incx;
                            x[ix - 1] = x[ix - 1] - temp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                    jx = jx - incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    if (x[j - 1] != zero) {
                        if (nounit) {
                            x[j - 1] = x[j - 1] / a[(j - 1) + (j - 1) * lda];
                        }
                        temp = x[j - 1];
                        for (i = j + 1; i <= n; i = i + 1) {
                            x[i - 1] = x[i - 1] - temp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    if (x[jx - 1] != zero) {
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] / a[(j - 1) + (j - 1) * lda];
                        }
                        temp = x[jx - 1];
                        ix = jx;
                        for (i = j + 1; i <= n; i = i + 1) {
                            ix += incx;
                            x[ix - 1] = x[ix - 1] - temp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                    jx += incx;
                }
            }
        }
    } else {
        //
        //        Form  x := inv( A**T )*x.
        //
        if (Mlsame_gmp(uplo, "U")) {
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[j - 1];
                    for (i = 1; i <= j - 1; i = i + 1) {
                        temp = temp - a[(i - 1) + (j - 1) * lda] * x[i - 1];
                    }
                    if (nounit) {
                        temp = temp / a[(j - 1) + (j - 1) * lda];
                    }
                    x[j - 1] = temp;
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[jx - 1];
                    ix = kx;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        temp = temp - a[(i - 1) + (j - 1) * lda] * x[ix - 1];
                        ix += incx;
                    }
                    if (nounit) {
                        temp = temp / a[(j - 1) + (j - 1) * lda];
                    }
                    x[jx - 1] = temp;
                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[j - 1];
                    for (i = n; i >= j + 1; i = i - 1) {
                        temp = temp - a[(i - 1) + (j - 1) * lda] * x[i - 1];
                    }
                    if (nounit) {
                        temp = temp / a[(j - 1) + (j - 1) * lda];
                    }
                    x[j - 1] = temp;
                }
            } else {
                kx += (n - 1) * incx;
                jx = kx;
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[jx - 1];
                    ix = kx;
                    for (i = n; i >= j + 1; i = i - 1) {
                        temp = temp - a[(i - 1) + (j - 1) * lda] * x[ix - 1];
                        ix = ix - incx;
                    }
                    if (nounit) {
                        temp = temp / a[(j - 1) + (j - 1) * lda];
                    }
                    x[jx - 1] = temp;
                    jx = jx - incx;
                }
            }
        }
    }
    //
    //     End of Rtrsv .
    //
}
