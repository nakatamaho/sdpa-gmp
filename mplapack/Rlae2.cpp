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

void Rlae2(dd_real const a, dd_real const b, dd_real const c, dd_real &rt1, dd_real &rt2) {
    //
    //     Compute the eigenvalues
    //
    dd_real sm = a + c;
    dd_real df = a - c;
    dd_real adf = abs(df);
    dd_real tb = b + b;
    dd_real ab = abs(tb);
    dd_real acmx = 0.0;
    dd_real acmn = 0.0;
    if (abs(a) > abs(c)) {
        acmx = a;
        acmn = c;
    } else {
        acmx = c;
        acmn = a;
    }
    const dd_real one = 1.0;
    dd_real rt = 0.0;
    const dd_real two = 2.0;
    if (adf > ab) {
        rt = adf * sqrt(one + pow2((ab / adf)));
    } else if (adf < ab) {
        rt = ab * sqrt(one + pow2((adf / ab)));
    } else {
        //
        //        Includes case AB=ADF=0
        //
        rt = ab * sqrt(two);
    }
    const dd_real zero = 0.0;
    const dd_real half = 0.5e0;
    if (sm < zero) {
        rt1 = half * (sm - rt);
        //
        //        Order of execution important.
        //        To get fully accurate smaller eigenvalue,
        //        next line needs to be executed in higher precision.
        //
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else if (sm > zero) {
        rt1 = half * (sm + rt);
        //
        //        Order of execution important.
        //        To get fully accurate smaller eigenvalue,
        //        next line needs to be executed in higher precision.
        //
        rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
    } else {
        //
        //        Includes case RT1 = RT2 = 0
        //
        rt1 = half * rt;
        rt2 = -half * rt;
    }
    //
    //     End of Rlae2
    //
}
