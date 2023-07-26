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

void Rlartg(dd_real const f, dd_real const g, dd_real &cs, dd_real &sn, dd_real &r) {
    dd_real safmin = 0.0;
    dd_real eps = 0.0;
    const dd_real two = 2.0;
    dd_real safmn2 = 0.0;
    const dd_real one = 1.0;
    dd_real safmx2 = 0.0;
    const dd_real zero = 0.0;
    dd_real f1 = 0.0;
    dd_real g1 = 0.0;
    dd_real scale = 0.0;
    mplapackint count = 0;
    mplapackint i = 0;
    //
    safmin = Rlamch_dd("S");
    eps = Rlamch_dd("E");
    safmn2 = pow(Rlamch_dd("B"), castINTEGER_dd(log(safmin / eps) / log(Rlamch_dd("B")) / two));
    safmx2 = one / safmn2;
    //        FIRST = .FALSE.
    //     END IF
    if (g == zero) {
        cs = one;
        sn = zero;
        r = f;
    } else if (f == zero) {
        cs = zero;
        sn = one;
        r = g;
    } else {
        f1 = f;
        g1 = g;
        scale = std::max(abs(f1), abs(g1));
        if (scale >= safmx2) {
            count = 0;
        statement_10:
            count++;
            f1 = f1 * safmn2;
            g1 = g1 * safmn2;
            scale = std::max(abs(f1), abs(g1));
            if (scale >= safmx2 && count < 20) {
                goto statement_10;
            }
            r = sqrt(pow2(f1) + pow2(g1));
            cs = f1 / r;
            sn = g1 / r;
            for (i = 1; i <= count; i = i + 1) {
                r = r * safmx2;
            }
        } else if (scale <= safmn2) {
            count = 0;
        statement_30:
            count++;
            f1 = f1 * safmx2;
            g1 = g1 * safmx2;
            scale = std::max(abs(f1), abs(g1));
            if (scale <= safmn2) {
                goto statement_30;
            }
            r = sqrt(pow2(f1) + pow2(g1));
            cs = f1 / r;
            sn = g1 / r;
            for (i = 1; i <= count; i = i + 1) {
                r = r * safmn2;
            }
        } else {
            r = sqrt(pow2(f1) + pow2(g1));
            cs = f1 / r;
            sn = g1 / r;
        }
        if (abs(f) > abs(g) && cs < zero) {
            cs = -cs;
            sn = -sn;
            r = -r;
        }
    }
    //
    //     End of Rlartg
    //
}
