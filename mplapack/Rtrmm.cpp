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

void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb) {
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
    bool lside = Mlsame_dd(side, "L");
    mplapackint nrowa = 0;
    if (lside) {
        nrowa = m;
    } else {
        nrowa = n;
    }
    bool nounit = Mlsame_dd(diag, "N");
    bool upper = Mlsame_dd(uplo, "U");
    //
    mplapackint info = 0;
    if ((!lside) && (!Mlsame_dd(side, "R"))) {
        info = 1;
    } else if ((!upper) && (!Mlsame_dd(uplo, "L"))) {
        info = 2;
    } else if ((!Mlsame_dd(transa, "N")) && (!Mlsame_dd(transa, "T")) && (!Mlsame_dd(transa, "C"))) {
        info = 3;
    } else if ((!Mlsame_dd(diag, "U")) && (!Mlsame_dd(diag, "N"))) {
        info = 4;
    } else if (m < 0) {
        info = 5;
    } else if (n < 0) {
        info = 6;
    } else if (lda < std::max((mplapackint)1, nrowa)) {
        info = 9;
    } else if (ldb < std::max((mplapackint)1, m)) {
        info = 11;
    }
    if (info != 0) {
        Mxerbla_dd("Rtrmm ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     And when  alpha.eq.zero.
    //
    const mpf_class zero = 0.0;
    mplapackint j = 0;
    mplapackint i = 0;
    if (alpha == zero) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = zero;
            }
        }
        return;
    }
    //
    //     Start the operations.
    //
    mplapackint k = 0;
    mpf_class temp = 0.0;
    const mpf_class one = 1.0;
    if (lside) {
        if (Mlsame_dd(transa, "N")) {
            //
            //           Form  B := alpha*A*B.
            //
            if (upper) {
                for (j = 1; j <= n; j = j + 1) {
                    for (k = 1; k <= m; k = k + 1) {
                        if (b[(k - 1) + (j - 1) * ldb] != zero) {
                            temp = alpha * b[(k - 1) + (j - 1) * ldb];
                            for (i = 1; i <= k - 1; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * a[(i - 1) + (k - 1) * lda];
                            }
                            if (nounit) {
                                temp = temp * a[(k - 1) + (k - 1) * lda];
                            }
                            b[(k - 1) + (j - 1) * ldb] = temp;
                        }
                    }
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (k = m; k >= 1; k = k - 1) {
                        if (b[(k - 1) + (j - 1) * ldb] != zero) {
                            temp = alpha * b[(k - 1) + (j - 1) * ldb];
                            b[(k - 1) + (j - 1) * ldb] = temp;
                            if (nounit) {
                                b[(k - 1) + (j - 1) * ldb] = b[(k - 1) + (j - 1) * ldb] * a[(k - 1) + (k - 1) * lda];
                            }
                            for (i = k + 1; i <= m; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * a[(i - 1) + (k - 1) * lda];
                            }
                        }
                    }
                }
            }
        } else {
            //
            //           Form  B := alpha*A**T*B.
            //
            if (upper) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = m; i >= 1; i = i - 1) {
                        temp = b[(i - 1) + (j - 1) * ldb];
                        if (nounit) {
                            temp = temp * a[(i - 1) + (i - 1) * lda];
                        }
                        for (k = 1; k <= i - 1; k = k + 1) {
                            temp += a[(k - 1) + (i - 1) * lda] * b[(k - 1) + (j - 1) * ldb];
                        }
                        b[(i - 1) + (j - 1) * ldb] = alpha * temp;
                    }
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= m; i = i + 1) {
                        temp = b[(i - 1) + (j - 1) * ldb];
                        if (nounit) {
                            temp = temp * a[(i - 1) + (i - 1) * lda];
                        }
                        for (k = i + 1; k <= m; k = k + 1) {
                            temp += a[(k - 1) + (i - 1) * lda] * b[(k - 1) + (j - 1) * ldb];
                        }
                        b[(i - 1) + (j - 1) * ldb] = alpha * temp;
                    }
                }
            }
        }
    } else {
        if (Mlsame_dd(transa, "N")) {
            //
            //           Form  B := alpha*B*A.
            //
            if (upper) {
                for (j = n; j >= 1; j = j - 1) {
                    temp = alpha;
                    if (nounit) {
                        temp = temp * a[(j - 1) + (j - 1) * lda];
                    }
                    for (i = 1; i <= m; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] = temp * b[(i - 1) + (j - 1) * ldb];
                    }
                    for (k = 1; k <= j - 1; k = k + 1) {
                        if (a[(k - 1) + (j - 1) * lda] != zero) {
                            temp = alpha * a[(k - 1) + (j - 1) * lda];
                            for (i = 1; i <= m; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * b[(i - 1) + (k - 1) * ldb];
                            }
                        }
                    }
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    temp = alpha;
                    if (nounit) {
                        temp = temp * a[(j - 1) + (j - 1) * lda];
                    }
                    for (i = 1; i <= m; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] = temp * b[(i - 1) + (j - 1) * ldb];
                    }
                    for (k = j + 1; k <= n; k = k + 1) {
                        if (a[(k - 1) + (j - 1) * lda] != zero) {
                            temp = alpha * a[(k - 1) + (j - 1) * lda];
                            for (i = 1; i <= m; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * b[(i - 1) + (k - 1) * ldb];
                            }
                        }
                    }
                }
            }
        } else {
            //
            //           Form  B := alpha*B*A**T.
            //
            if (upper) {
                for (k = 1; k <= n; k = k + 1) {
                    for (j = 1; j <= k - 1; j = j + 1) {
                        if (a[(j - 1) + (k - 1) * lda] != zero) {
                            temp = alpha * a[(j - 1) + (k - 1) * lda];
                            for (i = 1; i <= m; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * b[(i - 1) + (k - 1) * ldb];
                            }
                        }
                    }
                    temp = alpha;
                    if (nounit) {
                        temp = temp * a[(k - 1) + (k - 1) * lda];
                    }
                    if (temp != one) {
                        for (i = 1; i <= m; i = i + 1) {
                            b[(i - 1) + (k - 1) * ldb] = temp * b[(i - 1) + (k - 1) * ldb];
                        }
                    }
                }
            } else {
                for (k = n; k >= 1; k = k - 1) {
                    for (j = k + 1; j <= n; j = j + 1) {
                        if (a[(j - 1) + (k - 1) * lda] != zero) {
                            temp = alpha * a[(j - 1) + (k - 1) * lda];
                            for (i = 1; i <= m; i = i + 1) {
                                b[(i - 1) + (j - 1) * ldb] += temp * b[(i - 1) + (k - 1) * ldb];
                            }
                        }
                    }
                    temp = alpha;
                    if (nounit) {
                        temp = temp * a[(k - 1) + (k - 1) * lda];
                    }
                    if (temp != one) {
                        for (i = 1; i <= m; i = i + 1) {
                            b[(i - 1) + (k - 1) * ldb] = temp * b[(i - 1) + (k - 1) * ldb];
                        }
                    }
                }
            }
        }
    }
    //
    //     End of Rtrmm .
    //
}
