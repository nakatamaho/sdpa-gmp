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

void Rorgtr(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info) {
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
    //     Test the input arguments
    //
    info = 0;
    bool lquery = (lwork == -1);
    bool upper = Mlsame_dd(uplo, "U");
    if (!upper && !Mlsame_dd(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -4;
    } else if (lwork < std::max((mplapackint)1, n - 1) && !lquery) {
        info = -7;
    }
    //
    mplapackint nb = 0;
    mplapackint lwkopt = 0;
    if (info == 0) {
        if (upper) {
            nb = iMlaenv_dd(1, "Rorgql", " ", n - 1, n - 1, n - 1, -1);
        } else {
            nb = iMlaenv_dd(1, "Rorgqr", " ", n - 1, n - 1, n - 1, -1);
        }
        lwkopt = std::max((mplapackint)1, n - 1) * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla_dd("Rorgtr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    mplapackint j = 0;
    mplapackint i = 0;
    const mpf_class zero = 0.0;
    const mpf_class one = 1.0;
    mplapackint iinfo = 0;
    if (upper) {
        //
        //        Q was determined by a call to Rsytrd with UPLO = 'U'
        //
        //        Shift the vectors which define the elementary reflectors one
        //        column to the left, and set the last row and column of Q to
        //        those of the unit matrix
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j + 1) - 1) * lda];
            }
            a[(n - 1) + (j - 1) * lda] = zero;
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            a[(i - 1) + (n - 1) * lda] = zero;
        }
        a[(n - 1) + (n - 1) * lda] = one;
        //
        //        Generate Q(1:n-1,1:n-1)
        //
        Rorgql(n - 1, n - 1, n - 1, a, lda, tau, work, lwork, iinfo);
        //
    } else {
        //
        //        Q was determined by a call to Rsytrd with UPLO = 'L'.
        //
        //        Shift the vectors which define the elementary reflectors one
        //        column to the right, and set the first row and column of Q to
        //        those of the unit matrix
        //
        for (j = n; j >= 2; j = j - 1) {
            a[(j - 1) * lda] = zero;
            for (i = j + 1; i <= n; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j - 1) - 1) * lda];
            }
        }
        a[(1 - 1)] = one;
        for (i = 2; i <= n; i = i + 1) {
            a[(i - 1)] = zero;
        }
        if (n > 1) {
            //
            //           Generate Q(2:n,2:n)
            //
            Rorgqr(n - 1, n - 1, n - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, lwork, iinfo);
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Rorgtr
    //
}
