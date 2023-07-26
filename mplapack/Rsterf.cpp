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

void Rsterf(mplapackint const n, dd_real *d, dd_real *e, mplapackint &info) {
    dd_real eps = 0.0;
    dd_real eps2 = 0.0;
    dd_real safmin = 0.0;
    const dd_real one = 1.0;
    dd_real safmax = 0.0;
    const dd_real three = 3.0;
    dd_real ssfmax = 0.0;
    dd_real ssfmin = 0.0;
    dd_real rmax = 0.0;
    const mplapackint maxit = 30;
    mplapackint nmaxit = 0;
    const dd_real zero = 0.0;
    dd_real sigma = 0.0;
    mplapackint jtot = 0;
    mplapackint l1 = 0;
    mplapackint m = 0;
    mplapackint l = 0;
    mplapackint lsv = 0;
    mplapackint lend = 0;
    mplapackint lendsv = 0;
    dd_real anorm = 0.0;
    mplapackint iscale = 0;
    mplapackint i = 0;
    dd_real p = 0.0;
    dd_real rte = 0.0;
    dd_real rt1 = 0.0;
    dd_real rt2 = 0.0;
    const dd_real two = 2.0;
    dd_real r = 0.0;
    dd_real c = 0.0;
    dd_real s = 0.0;
    dd_real gamma = 0.0;
    dd_real bb = 0.0;
    dd_real oldc = 0.0;
    dd_real oldgam = 0.0;
    dd_real alpha = 0.0;
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
    //     Test the input parameters.
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n < 0) {
        info = -1;
        Mxerbla_dd("Rsterf", -info);
        return;
    }
    if (n <= 1) {
        return;
    }
    //
    //     Determine the unit roundoff for this environment.
    //
    eps = Rlamch_dd("E");
    eps2 = pow2(eps);
    safmin = Rlamch_dd("S");
    safmax = one / safmin;
    ssfmax = sqrt(safmax) / three;
    ssfmin = sqrt(safmin) / eps2;
    rmax = Rlamch_dd("O");
    //
    //     Compute the eigenvalues of the tridiagonal matrix.
    //
    nmaxit = n * maxit;
    sigma = zero;
    jtot = 0;
    //
    //     Determine where the matrix splits and choose QL or QR iteration
    //     for each block, according to whether top or bottom diagonal
    //     element is smaller.
    //
    l1 = 1;
//
statement_10:
    if (l1 > n) {
        goto statement_170;
    }
    if (l1 > 1) {
        e[(l1 - 1) - 1] = zero;
    }
    for (m = l1; m <= n - 1; m = m + 1) {
        if (abs(e[m - 1]) <= (sqrt(abs(d[m - 1])) * sqrt(abs(d[(m + 1) - 1]))) * eps) {
            e[m - 1] = zero;
            goto statement_30;
        }
    }
    m = n;
//
statement_30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
        goto statement_10;
    }
    //
    //     Scale submatrix in rows and columns L to LEND
    //
    anorm = Rlanst("M", lend - l + 1, &d[l - 1], &e[l - 1]);
    iscale = 0;
    if (anorm == zero) {
        goto statement_10;
    }
    if ((anorm > ssfmax)) {
        iscale = 1;
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l + 1, 1, &d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmax, lend - l, 1, &e[l - 1], n, info);
    } else if (anorm < ssfmin) {
        iscale = 2;
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l + 1, 1, &d[l - 1], n, info);
        Rlascl("G", 0, 0, anorm, ssfmin, lend - l, 1, &e[l - 1], n, info);
    }
    //
    for (i = l; i <= lend - 1; i = i + 1) {
        e[i - 1] = pow2(e[i - 1]);
    }
    //
    //     Choose between QL and QR iteration
    //
    if (abs(d[lend - 1]) < abs(d[l - 1])) {
        lend = lsv;
        l = lendsv;
    }
    //
    if (lend >= l) {
    //
    //        QL Iteration
    //
    //        Look for small subdiagonal element.
    //
    statement_50:
        if (l != lend) {
            for (m = l; m <= lend - 1; m = m + 1) {
                if (abs(e[m - 1]) <= eps2 * abs(d[m - 1] * d[(m + 1) - 1])) {
                    goto statement_70;
                }
            }
        }
        m = lend;
    //
    statement_70:
        if (m < lend) {
            e[m - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_90;
        }
        //
        //        If remaining matrix is 2 by 2, use Rlae2 to compute its
        //        eigenvalues.
        //
        if (m == l + 1) {
            rte = sqrt(e[l - 1]);
            Rlae2(d[l - 1], rte, d[(l + 1) - 1], rt1, rt2);
            d[l - 1] = rt1;
            d[(l + 1) - 1] = rt2;
            e[l - 1] = zero;
            l += 2;
            if (l <= lend) {
                goto statement_50;
            }
            goto statement_150;
        }
        //
        if (jtot == nmaxit) {
            goto statement_150;
        }
        jtot++;
        //
        //        Form shift.
        //
        rte = sqrt(e[l - 1]);
        sigma = (d[(l + 1) - 1] - p) / (two * rte);
        r = Rlapy2(sigma, one);
        sigma = p - (rte / (sigma + sign(r, sigma)));
        //
        c = one;
        s = zero;
        gamma = d[m - 1] - sigma;
        p = gamma * gamma;
        //
        //        Inner loop
        //
        for (i = m - 1; i >= l; i = i - 1) {
            bb = e[i - 1];
            r = p + bb;
            if (i != m - 1) {
                e[(i + 1) - 1] = s * r;
            }
            oldc = c;
            c = p / r;
            s = bb / r;
            oldgam = gamma;
            alpha = d[i - 1];
            gamma = c * (alpha - sigma) - s * oldgam;
            d[(i + 1) - 1] = oldgam + (alpha - gamma);
            if (c != zero) {
                p = (gamma * gamma) / c;
            } else {
                p = oldc * bb;
            }
        }
        //
        e[l - 1] = s * p;
        d[l - 1] = sigma + gamma;
        goto statement_50;
    //
    //        Eigenvalue found.
    //
    statement_90:
        d[l - 1] = p;
        //
        l++;
        if (l <= lend) {
            goto statement_50;
        }
        goto statement_150;
        //
    } else {
    //
    //        QR Iteration
    //
    //        Look for small superdiagonal element.
    //
    statement_100:
        for (m = l; m >= lend + 1; m = m - 1) {
            if (abs(e[(m - 1) - 1]) <= eps2 * abs(d[m - 1] * d[(m - 1) - 1])) {
                goto statement_120;
            }
        }
        m = lend;
    //
    statement_120:
        if (m > lend) {
            e[(m - 1) - 1] = zero;
        }
        p = d[l - 1];
        if (m == l) {
            goto statement_140;
        }
        //
        //        If remaining matrix is 2 by 2, use Rlae2 to compute its
        //        eigenvalues.
        //
        if (m == l - 1) {
            rte = sqrt(e[(l - 1) - 1]);
            Rlae2(d[l - 1], rte, d[(l - 1) - 1], rt1, rt2);
            d[l - 1] = rt1;
            d[(l - 1) - 1] = rt2;
            e[(l - 1) - 1] = zero;
            l = l - 2;
            if (l >= lend) {
                goto statement_100;
            }
            goto statement_150;
        }
        //
        if (jtot == nmaxit) {
            goto statement_150;
        }
        jtot++;
        //
        //        Form shift.
        //
        rte = sqrt(e[(l - 1) - 1]);
        sigma = (d[(l - 1) - 1] - p) / (two * rte);
        r = Rlapy2(sigma, one);
        sigma = p - (rte / (sigma + sign(r, sigma)));
        //
        c = one;
        s = zero;
        gamma = d[m - 1] - sigma;
        p = gamma * gamma;
        //
        //        Inner loop
        //
        for (i = m; i <= l - 1; i = i + 1) {
            bb = e[i - 1];
            r = p + bb;
            if (i != m) {
                e[(i - 1) - 1] = s * r;
            }
            oldc = c;
            c = p / r;
            s = bb / r;
            oldgam = gamma;
            alpha = d[(i + 1) - 1];
            gamma = c * (alpha - sigma) - s * oldgam;
            d[i - 1] = oldgam + (alpha - gamma);
            if (c != zero) {
                p = (gamma * gamma) / c;
            } else {
                p = oldc * bb;
            }
        }
        //
        e[(l - 1) - 1] = s * p;
        d[l - 1] = sigma + gamma;
        goto statement_100;
    //
    //        Eigenvalue found.
    //
    statement_140:
        d[l - 1] = p;
        //
        l = l - 1;
        if (l >= lend) {
            goto statement_100;
        }
        goto statement_150;
        //
    }
//
//     Undo scaling if necessary
//
statement_150:
    if (iscale == 1) {
        Rlascl("G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n, info);
    }
    if (iscale == 2) {
        Rlascl("G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, &d[lsv - 1], n, info);
    }
    //
    //     Check for no convergence to an eigenvalue after a total
    //     of N*MAXIT iterations.
    //
    if (jtot < nmaxit) {
        goto statement_10;
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        if (e[i - 1] != zero) {
            info++;
        }
    }
    goto statement_180;
//
//     Sort eigenvalues in increasing order.
//
statement_170:
    Rlasrt("I", n, d, info);
//
statement_180:;
    //
    //     End of Rsterf
    //
}
