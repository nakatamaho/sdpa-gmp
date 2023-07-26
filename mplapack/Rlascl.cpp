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

void Rlascl(const char *type, mplapackint const kl, mplapackint const ku, mpf_class const cfrom, mpf_class const cto, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info) {
    mplapackint itype = 0;
    const mpf_class zero = 0.0;
    mpf_class smlnum = 0.0;
    const mpf_class one = 1.0;
    mpf_class bignum = 0.0;
    mpf_class cfromc = 0.0;
    mpf_class ctoc = 0.0;
    mpf_class cfrom1 = 0.0;
    mpf_class mul = 0.0;
    bool done = false;
    mpf_class cto1 = 0.0;
    mplapackint j = 0;
    mplapackint i = 0;
    mplapackint k3 = 0;
    mplapackint k4 = 0;
    mplapackint k1 = 0;
    mplapackint k2 = 0;
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    //
    if (Mlsame_dd(type, "G")) {
        itype = 0;
    } else if (Mlsame_dd(type, "L")) {
        itype = 1;
    } else if (Mlsame_dd(type, "U")) {
        itype = 2;
    } else if (Mlsame_dd(type, "H")) {
        itype = 3;
    } else if (Mlsame_dd(type, "B")) {
        itype = 4;
    } else if (Mlsame_dd(type, "Q")) {
        itype = 5;
    } else if (Mlsame_dd(type, "Z")) {
        itype = 6;
    } else {
        itype = -1;
    }
    //
    if (itype == -1) {
        info = -1;
    } else if (cfrom == zero || Risnan(cfrom)) {
        info = -4;
    } else if (Risnan(cto)) {
        info = -5;
    } else if (m < 0) {
        info = -6;
    } else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
        info = -7;
    } else if (itype <= 3 && lda < std::max((mplapackint)1, m)) {
        info = -9;
    } else if (itype >= 4) {
        if (kl < 0 || kl > std::max(m - 1, (mplapackint)0)) {
            info = -2;
        } else if (ku < 0 || ku > std::max(n - 1, (mplapackint)0) || ((itype == 4 || itype == 5) && kl != ku)) {
            info = -3;
        } else if ((itype == 4 && lda < kl + 1) || (itype == 5 && lda < ku + 1) || (itype == 6 && lda < 2 * kl + ku + 1)) {
            info = -9;
        }
    }
    //
    if (info != 0) {
        Mxerbla_dd("Rlascl", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || m == 0) {
        return;
    }
    //
    //     Get machine parameters
    //
    smlnum = Rlamch_dd("S");
    bignum = one / smlnum;
    //
    cfromc = cfrom;
    ctoc = cto;
//
statement_10:
    cfrom1 = cfromc * smlnum;
    if (cfrom1 == cfromc) {
        //        CFROMC is an inf.  Multiply by a correctly signed zero for
        //        finite CTOC, or a NaN if CTOC is infinite.
        mul = ctoc / cfromc;
        done = true;
        cto1 = ctoc;
    } else {
        cto1 = ctoc / bignum;
        if (cto1 == ctoc) {
            //           CTOC is either 0 or an inf.  In both cases, CTOC itself
            //           serves as the correct multiplication factor.
            mul = ctoc;
            done = true;
            cfromc = one;
        } else if (abs(cfrom1) > abs(ctoc) && ctoc != zero) {
            mul = smlnum;
            done = false;
            cfromc = cfrom1;
        } else if (abs(cto1) > abs(cfromc)) {
            mul = bignum;
            done = false;
            ctoc = cto1;
        } else {
            mul = ctoc / cfromc;
            done = true;
        }
    }
    //
    if (itype == 0) {
        //
        //        Full matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 1) {
        //
        //        Lower triangular matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 2) {
        //
        //        Upper triangular matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= std::min(j, m); i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 3) {
        //
        //        Upper Hessenberg matrix
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= std::min(j + 1, m); i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 4) {
        //
        //        Lower half of a symmetric band matrix
        //
        k3 = kl + 1;
        k4 = n + 1;
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= std::min(k3, k4 - j); i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 5) {
        //
        //        Upper half of a symmetric band matrix
        //
        k1 = ku + 2;
        k3 = ku + 1;
        for (j = 1; j <= n; j = j + 1) {
            for (i = std::max(k1 - j, (mplapackint)1); i <= k3; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    } else if (itype == 6) {
        //
        //        Band matrix
        //
        k1 = kl + ku + 2;
        k2 = kl + 1;
        k3 = 2 * kl + ku + 1;
        k4 = kl + ku + 1 + m;
        for (j = 1; j <= n; j = j + 1) {
            for (i = std::max(k1 - j, k2); i <= std::min(k3, k4 - j); i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * mul;
            }
        }
        //
    }
    //
    if (!done) {
        goto statement_10;
    }
    //
    //     End of Rlascl
    //
}
