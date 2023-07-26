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

void Rsytrd(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *d, mpf_class *e, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info) {
    //
    //     Test the input parameters
    //
    info = 0;
    bool upper = Mlsame_dd(uplo, "U");
    bool lquery = (lwork == -1);
    if (!upper && !Mlsame_dd(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -4;
    } else if (lwork < 1 && !lquery) {
        info = -9;
    }
    //
    mplapackint nb = 0;
    mplapackint lwkopt = 0;
    if (info == 0) {
        //
        //        Determine the block size.
        //
        nb = iMlaenv_dd(1, "Rsytrd", uplo, n, -1, -1, -1);
        lwkopt = n * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla_dd("Rsytrd", -info);
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
    mplapackint nx = n;
    mplapackint iws = 1;
    mplapackint ldwork = 0;
    mplapackint nbmin = 0;
    if (nb > 1 && nb < n) {
        //
        //        Determine when to cross over from blocked to unblocked code
        //        (last block is always handled by unblocked code).
        //
        nx = std::max(nb, iMlaenv_dd(3, "Rsytrd", uplo, n, -1, -1, -1));
        if (nx < n) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            ldwork = n;
            iws = ldwork * nb;
            if (lwork < iws) {
                //
                //              Not enough workspace to use optimal NB:  determine the
                //              minimum value of NB, and reduce NB or force use of
                //              unblocked code by setting NX = N.
                //
                nb = std::max(lwork / ldwork, (mplapackint)1);
                nbmin = iMlaenv_dd(2, "Rsytrd", uplo, n, -1, -1, -1);
                if (nb < nbmin) {
                    nx = n;
                }
            }
        } else {
            nx = n;
        }
    } else {
        nb = 1;
    }
    //
    mplapackint kk = 0;
    mplapackint i = 0;
    const mpf_class one = 1.0;
    mplapackint j = 0;
    mplapackint iinfo = 0;
    if (upper) {
        //
        //        Reduce the upper triangle of A.
        //        Columns 1:kk are handled by the unblocked method.
        //
        kk = n - ((n - nx + nb - 1) / nb) * nb;
        for (i = n - nb + 1; i >= kk + 1; i = i - nb) {
            //
            //           Reduce columns i:i+nb-1 to tridiagonal form and form the
            //           matrix W which is needed to update the unreduced part of
            //           the matrix
            //
            Rlatrd(uplo, i + nb - 1, nb, a, lda, e, tau, work, ldwork);
            //
            //           Update the unreduced submatrix A(1:i-1,1:i-1), using an
            //           update of the form:  A := A - V*W**T - W*V**T
            //
            Rsyr2k(uplo, "No transpose", i - 1, nb, -one, &a[(i - 1) * lda], lda, work, ldwork, one, a, lda);
            //
            //           Copy superdiagonal elements back into A, and diagonal
            //           elements into D
            //
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[((j - 1) - 1) + (j - 1) * lda] = e[(j - 1) - 1];
                d[j - 1] = a[(j - 1) + (j - 1) * lda];
            }
        }
        //
        //        Use unblocked code to reduce the last or only block
        //
        Rsytd2(uplo, kk, a, lda, d, e, tau, iinfo);
    } else {
        //
        //        Reduce the lower triangle of A
        //
        for (i = 1; i <= n - nx; i = i + nb) {
            //
            //           Reduce columns i:i+nb-1 to tridiagonal form and form the
            //           matrix W which is needed to update the unreduced part of
            //           the matrix
            //
            Rlatrd(uplo, n - i + 1, nb, &a[(i - 1) + (i - 1) * lda], lda, &e[i - 1], &tau[i - 1], work, ldwork);
            //
            //           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
            //           an update of the form:  A := A - V*W**T - W*V**T
            //
            Rsyr2k(uplo, "No transpose", n - i - nb + 1, nb, -one, &a[((i + nb) - 1) + (i - 1) * lda], lda, &work[(nb + 1) - 1], ldwork, one, &a[((i + nb) - 1) + ((i + nb) - 1) * lda], lda);
            //
            //           Copy subdiagonal elements back into A, and diagonal
            //           elements into D
            //
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[((j + 1) - 1) + (j - 1) * lda] = e[j - 1];
                d[j - 1] = a[(j - 1) + (j - 1) * lda];
            }
        }
        //
        //        Use unblocked code to reduce the last or only block
        //
        Rsytd2(uplo, n - i + 1, &a[(i - 1) + (i - 1) * lda], lda, &d[i - 1], &e[i - 1], &tau[i - 1], iinfo);
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rsytrd
    //
}
