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

void Rlasrt(const char *id, mplapackint const n, mpf_class *d, mplapackint &info) {
    mplapackint dir = 0;
    mplapackint stkpnt = 0;
    mplapackint stacklen = 32;
    mplapackint stack[2 * stacklen];
    mplapackint ldstack = 2;
    mplapackint start = 0;
    mplapackint endd = 0;
    const mplapackint select = 20;
    mplapackint i = 0;
    mplapackint j = 0;
    mpf_class dmnmx = 0.0;
    mpf_class d1 = 0.0;
    mpf_class d2 = 0.0;
    mpf_class d3 = 0.0;
    mpf_class tmp = 0.0;
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    dir = -1;
    if (Mlsame_dd(id, "D")) {
        dir = 0;
    } else if (Mlsame_dd(id, "I")) {
        dir = 1;
    }
    if (dir == -1) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    }
    if (info != 0) {
        Mxerbla_dd("Rlasrt", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    stkpnt = 1;
    stack[(1 - 1) + (1 - 1) * ldstack] = 1;
    stack[(2 - 1) + (1 - 1) * ldstack] = n;
statement_10:
    start = stack[(1 - 1) + (stkpnt - 1) * ldstack];
    endd = stack[(2 - 1) + (stkpnt - 1) * ldstack];
    stkpnt = stkpnt - 1;
    if (endd - start <= select && endd - start > 0) {
        //
        //        Do Insertion sort on D( START:ENDD )
        //
        if (dir == 0) {
            //
            //           Sort into decreasing order
            //
            for (i = start + 1; i <= endd; i = i + 1) {
                for (j = i; j >= start + 1; j = j - 1) {
                    if (d[j - 1] > d[(j - 1) - 1]) {
                        dmnmx = d[j - 1];
                        d[j - 1] = d[(j - 1) - 1];
                        d[(j - 1) - 1] = dmnmx;
                    } else {
                        goto statement_30;
                    }
                }
            statement_30:;
            }
            //
        } else {
            //
            //           Sort into increasing order
            //
            for (i = start + 1; i <= endd; i = i + 1) {
                for (j = i; j >= start + 1; j = j - 1) {
                    if (d[j - 1] < d[(j - 1) - 1]) {
                        dmnmx = d[j - 1];
                        d[j - 1] = d[(j - 1) - 1];
                        d[(j - 1) - 1] = dmnmx;
                    } else {
                        goto statement_50;
                    }
                }
            statement_50:;
            }
            //
        }
        //
    } else if (endd - start > select) {
        //
        //        Partition D( START:ENDD ) and stack parts, largest one first
        //
        //        Choose partition entry as median of 3
        //
        d1 = d[start - 1];
        d2 = d[endd - 1];
        i = (start + endd) / 2;
        d3 = d[i - 1];
        if (d1 < d2) {
            if (d3 < d1) {
                dmnmx = d1;
            } else if (d3 < d2) {
                dmnmx = d3;
            } else {
                dmnmx = d2;
            }
        } else {
            if (d3 < d2) {
                dmnmx = d2;
            } else if (d3 < d1) {
                dmnmx = d3;
            } else {
                dmnmx = d1;
            }
        }
        //
        if (dir == 0) {
            //
            //           Sort into decreasing order
            //
            i = start - 1;
            j = endd + 1;
        statement_60:
        statement_70:
            j = j - 1;
            if (d[j - 1] < dmnmx) {
                goto statement_70;
            }
        statement_80:
            i++;
            if (d[i - 1] > dmnmx) {
                goto statement_80;
            }
            if (i < j) {
                tmp = d[i - 1];
                d[i - 1] = d[j - 1];
                d[j - 1] = tmp;
                goto statement_60;
            }
            if (j - start > endd - j - 1) {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = start;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = j;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = endd;
            } else {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = endd;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = start;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = j;
            }
        } else {
            //
            //           Sort into increasing order
            //
            i = start - 1;
            j = endd + 1;
        statement_90:
        statement_100:
            j = j - 1;
            if (d[j - 1] > dmnmx) {
                goto statement_100;
            }
        statement_110:
            i++;
            if (d[i - 1] < dmnmx) {
                goto statement_110;
            }
            if (i < j) {
                tmp = d[i - 1];
                d[i - 1] = d[j - 1];
                d[j - 1] = tmp;
                goto statement_90;
            }
            if (j - start > endd - j - 1) {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = start;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = j;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = endd;
            } else {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = endd;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * ldstack] = start;
                stack[(2 - 1) + (stkpnt - 1) * ldstack] = j;
            }
        }
    }
    if (stkpnt > 0) {
        goto statement_10;
    }
    //
    //     End of Rlasrt
    //
}
