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

void Rscal(mplapackint const n, mpf_class const da, mpf_class *dx, mplapackint const incx) {
    if (n <= 0 || incx <= 0) {
        return;
    }

    mplapackint m = 0;
    mplapackint i = 0;
    mplapackint mp1 = 0;
    mplapackint nincx = 0;

    if (incx == 1) {
        // Code for increment equal to 1
        m = n % 5; // Using modulo operator for cleanup
        if (m != 0) {
            for (i = 0; i < m; i++) {
                mpf_class temp = dx[i];
                temp *= da;
                dx[i] = temp;
            }
            if (n < 5) {
                return;
            }
        }
        for (i = m; i < n; i += 5) {
            mpf_class temp1 = dx[i];
            temp1 *= da;
            dx[i] = temp1;

            mpf_class temp2 = dx[i + 1];
            temp2 *= da;
            dx[i + 1] = temp2;

            mpf_class temp3 = dx[i + 2];
            temp3 *= da;
            dx[i + 2] = temp3;

            mpf_class temp4 = dx[i + 3];
            temp4 *= da;
            dx[i + 3] = temp4;

            mpf_class temp5 = dx[i + 4];
            temp5 *= da;
            dx[i + 4] = temp5;
        }
    } else {
        // Code for increment not equal to 1
        mplapackint ix = 0;
        for (i = 0; i < n; i++) {
            mpf_class temp = dx[ix];
            temp *= da;
            dx[ix] = temp;
            ix += incx;
        }
    }
}
