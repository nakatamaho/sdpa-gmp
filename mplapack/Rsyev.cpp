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
#include <mplapack_gmp.h>

void Rsyev(const char *jobz, const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *w, mpf_class *work, mplapackint const lwork, mplapackint &info) {
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame_gmp(jobz, "V");
    bool lower = Mlsame_gmp(uplo, "L");
    bool lquery = (lwork == -1);
    //
    info = 0;
    if (!(wantz || Mlsame_gmp(jobz, "N"))) {
        info = -1;
    } else if (!(lower || Mlsame_gmp(uplo, "U"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < std::max((mplapackint)1, n)) {
        info = -5;
    }
    //
    mplapackint nb = 0;
    mplapackint lwkopt = 0;
    if (info == 0) {
        nb = iMlaenv_gmp(1, "Rsytrd", uplo, n, -1, -1, -1);
        lwkopt = std::max((mplapackint)1, (nb + 2) * n);
        work[1 - 1] = lwkopt;
        //
        if (lwork < std::max((mplapackint)1, 3 * n - 1) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla_gmp("Rsyev", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    const mpf_class one = 1.0;
    if (n == 1) {
        w[1 - 1] = a[(1 - 1)];
        work[1 - 1] = 2;
        if (wantz) {
            a[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    mpf_class safmin = Rlamch_gmp("Safe minimum");
    mpf_class eps = Rlamch_gmp("Precision");
    mpf_class smlnum = safmin / eps;
    mpf_class bignum = one / smlnum;
    mpf_class rmin = sqrt(smlnum);
    mpf_class rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    mpf_class anrm = Rlansy("M", uplo, n, a, lda, work);
    mplapackint iscale = 0;
    const mpf_class zero = 0.0;
    mpf_class sigma = 0.0;
    if (anrm > zero && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        Rlascl(uplo, 0, 0, one, sigma, n, n, a, lda, info);
    }
    //
    //     Call Rsytrd to reduce symmetric matrix to tridiagonal form.
    //
    mplapackint inde = 1;
    mplapackint indtau = inde + n;
    mplapackint indwrk = indtau + n;
    mplapackint llwork = lwork - indwrk + 1;
    mplapackint iinfo = 0;
    Rsytrd(uplo, n, a, lda, w, &work[inde - 1], &work[indtau - 1], &work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Rorgtr to generate the orthogonal matrix, then call Rsteqr.
    //
    if (!wantz) {
        Rsterf(n, w, &work[inde - 1], info);
    } else {
        Rorgtr(uplo, n, a, lda, &work[indtau - 1], &work[indwrk - 1], llwork, iinfo);
        Rsteqr(jobz, n, w, &work[inde - 1], a, lda, &work[indtau - 1], info);
    }
    //
    //     If matrix was scaled, then rescale eigenvalues appropriately.
    //
    mplapackint imax = 0;
    if (iscale == 1) {
        if (info == 0) {
            imax = n;
        } else {
            imax = info - 1;
        }
        Rscal(imax, one / sigma, w, 1);
    }
    //
    //     Set WORK(1) to optimal workspace size.
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rsyev
    //
}
