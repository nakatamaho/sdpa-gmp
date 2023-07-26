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

void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc) {
    //
    //  -- Reference BLAS level3 routine --
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Parameters ..
    //     ..
    //
    //     Test the input parameters.
    //
    mplapackint nrowa = 0;
    if (Mlsame_gmp(trans, "N")) {
        nrowa = n;
    } else {
        nrowa = k;
    }
    bool upper = Mlsame_gmp(uplo, "U");
    //
    mplapackint info = 0;
    if ((!upper) && (!Mlsame_gmp(uplo, "L"))) {
        info = 1;
    } else if ((!Mlsame_gmp(trans, "N")) && (!Mlsame_gmp(trans, "T")) && (!Mlsame_gmp(trans, "C"))) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (k < 0) {
        info = 4;
    } else if (lda < std::max((mplapackint)1, nrowa)) {
        info = 7;
    } else if (ldb < std::max((mplapackint)1, nrowa)) {
        info = 9;
    } else if (ldc < std::max((mplapackint)1, n)) {
        info = 12;
    }
    if (info != 0) {
        Mxerbla_gmp("Rsyr2k", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const mpf_class zero = 0.0;
    const mpf_class one = 1.0;
    if ((n == 0) || (((alpha == zero) || (k == 0)) && (beta == one))) {
        return;
    }
    //
    //     And when  alpha.eq.zero.
    //
    mplapackint j = 0;
    mplapackint i = 0;
    if (alpha == zero) {
        if (upper) {
            if (beta == zero) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        } else {
            if (beta == zero) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j; i <= n; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j; i <= n; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        }
        return;
    }
    //
    //     Start the operations.
    //
    mplapackint l = 0;
    mpf_class temp1 = 0.0;
    mpf_class temp2 = 0.0;
    if (Mlsame_gmp(trans, "N")) {
        //
        //        Form  C := alpha*A*B**T + alpha*B*A**T + C.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                if (beta == zero) {
                    for (i = 1; i <= j; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = 1; i <= j; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    if ((a[(j - 1) + (l - 1) * lda] != zero) || (b[(j - 1) + (l - 1) * ldb] != zero)) {
                        temp1 = alpha * b[(j - 1) + (l - 1) * ldb];
                        temp2 = alpha * a[(j - 1) + (l - 1) * lda];
                        for (i = 1; i <= j; i = i + 1) {
                            c[(i - 1) + (j - 1) * ldc] += a[(i - 1) + (l - 1) * lda] * temp1 + b[(i - 1) + (l - 1) * ldb] * temp2;
                        }
                    }
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                if (beta == zero) {
                    for (i = j; i <= n; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = j; i <= n; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    if ((a[(j - 1) + (l - 1) * lda] != zero) || (b[(j - 1) + (l - 1) * ldb] != zero)) {
                        temp1 = alpha * b[(j - 1) + (l - 1) * ldb];
                        temp2 = alpha * a[(j - 1) + (l - 1) * lda];
                        for (i = j; i <= n; i = i + 1) {
                            c[(i - 1) + (j - 1) * ldc] += a[(i - 1) + (l - 1) * lda] * temp1 + b[(i - 1) + (l - 1) * ldb] * temp2;
                        }
                    }
                }
            }
        }
    } else {
        //
        //        Form  C := alpha*A**T*B + alpha*B**T*A + C.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    temp1 = zero;
                    temp2 = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp1 += a[(l - 1) + (i - 1) * lda] * b[(l - 1) + (j - 1) * ldb];
                        temp2 += b[(l - 1) + (i - 1) * ldb] * a[(l - 1) + (j - 1) * lda];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp1 + alpha * temp2;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc] + alpha * temp1 + alpha * temp2;
                    }
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = j; i <= n; i = i + 1) {
                    temp1 = zero;
                    temp2 = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp1 += a[(l - 1) + (i - 1) * lda] * b[(l - 1) + (j - 1) * ldb];
                        temp2 += b[(l - 1) + (i - 1) * ldb] * a[(l - 1) + (j - 1) * lda];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp1 + alpha * temp2;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc] + alpha * temp1 + alpha * temp2;
                    }
                }
            }
        }
    }
    //
    //     End of Rsyr2k.
    //
}
