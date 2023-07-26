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

void Rorgql(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (k < 0 || k > n) {
        info = -3;
    } else if (lda < std::max((mplapackint)1, m)) {
        info = -5;
    }
    //
    mplapackint lwkopt = 0;
    mplapackint nb = 0;
    if (info == 0) {
        if (n == 0) {
            lwkopt = 1;
        } else {
            nb = iMlaenv_dd(1, "Rorgql", " ", m, n, k, -1);
            lwkopt = n * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < std::max((mplapackint)1, n) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla_dd("Rorgql", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    mplapackint nbmin = 2;
    mplapackint nx = 0;
    mplapackint iws = n;
    mplapackint ldwork = 0;
    if (nb > 1 && nb < k) {
        //
        //        Determine when to cross over from blocked to unblocked code.
        //
        nx = std::max((mplapackint)0, iMlaenv_dd(3, "Rorgql", " ", m, n, k, -1));
        if (nx < k) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            ldwork = n;
            iws = ldwork * nb;
            if (lwork < iws) {
                //
                //              Not enough workspace to use optimal NB:  reduce NB and
                //              determine the minimum value of NB.
                //
                nb = lwork / ldwork;
                nbmin = std::max((mplapackint)2, iMlaenv_dd(2, "Rorgql", " ", m, n, k, -1));
            }
        }
    }
    //
    mplapackint kk = 0;
    mplapackint j = 0;
    mplapackint i = 0;
    const mpf_class zero = 0.0;
    if (nb >= nbmin && nb < k && nx < k) {
        //
        //        Use blocked code after the first block.
        //        The last kk columns are handled by the block method.
        //
        kk = std::min(k, ((k - nx + nb - 1) / nb) * nb);
        //
        //        Set A(m-kk+1:m,1:n-kk) to zero.
        //
        for (j = 1; j <= n - kk; j = j + 1) {
            for (i = m - kk + 1; i <= m; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = zero;
            }
        }
    } else {
        kk = 0;
    }
    //
    //     Use unblocked code for the first or only block.
    //
    mplapackint iinfo = 0;
    Rorg2l(m - kk, n - kk, k - kk, a, lda, tau, work, iinfo);
    //
    mplapackint ib = 0;
    mplapackint l = 0;
    if (kk > 0) {
        //
        //        Use blocked code
        //
        for (i = k - kk + 1; i <= k; i = i + nb) {
            ib = std::min(nb, k - i + 1);
            if (n - k + i > 1) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i+ib-1) . . . H(i+1) H(i)
                //
                Rlarft("Backward", "Columnwise", m - k + i + ib - 1, ib, &a[((n - k + i) - 1) * lda], lda, &tau[i - 1], work, ldwork);
                //
                //              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
                //
                Rlarfb("Left", "No transpose", "Backward", "Columnwise", m - k + i + ib - 1, n - k + i - 1, ib, &a[((n - k + i) - 1) * lda], lda, work, ldwork, a, lda, &work[(ib + 1) - 1], ldwork);
            }
            //
            //           Apply H to rows 1:m-k+i+ib-1 of current block
            //
            Rorg2l(m - k + i + ib - 1, ib, ib, &a[((n - k + i) - 1) * lda], lda, &tau[i - 1], work, iinfo);
            //
            //           Set rows m-k+i+ib:m of current block to zero
            //
            for (j = n - k + i; j <= n - k + i + ib - 1; j = j + 1) {
                for (l = m - k + i + ib; l <= m; l = l + 1) {
                    a[(l - 1) + (j - 1) * lda] = zero;
                }
            }
        }
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Rorgql
    //
}
