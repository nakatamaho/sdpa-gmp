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

void Rlarf(const char *side, mplapackint const m, mplapackint const n, mpf_class *v, mplapackint const incv, mpf_class const tau, mpf_class *c, mplapackint const ldc, mpf_class *work) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    bool applyleft = Mlsame_dd(side, "L");
    mplapackint lastv = 0;
    mplapackint lastc = 0;
    const mpf_class zero = 0.0;
    mplapackint i = 0;
    if (tau != zero) {
        //     Set up variables for scanning V.  LASTV begins pointing to the end
        //     of V.
        if (applyleft) {
            lastv = m;
        } else {
            lastv = n;
        }
        if (incv > 0) {
            i = 1 + (lastv - 1) * incv;
        } else {
            i = 1;
        }
        //     Look for the last non-zero row in V.
        while (lastv > 0 && v[i - 1] == zero) {
            lastv = lastv - 1;
            i = i - incv;
        }
        if (applyleft) {
            //     Scan for the last non-zero column in C(1:lastv,:).
            lastc = iMladlc(lastv, n, c, ldc);
        } else {
            //     Scan for the last non-zero row in C(:,1:lastv).
            lastc = iMladlr(m, lastv, c, ldc);
        }
    }
    //     Note that lastc.eq.0 renders the BLAS operations null; no special
    //     case is needed at this level.
    const mpf_class one = 1.0;
    if (applyleft) {
        //
        //        Form  H * C
        //
        if (lastv > 0) {
            //
            //           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
            //
            Rgemv("Transpose", lastv, lastc, one, c, ldc, v, incv, zero, work, 1);
            //
            //           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
            //
            Rger(lastv, lastc, -tau, v, incv, work, 1, c, ldc);
        }
    } else {
        //
        //        Form  C * H
        //
        if (lastv > 0) {
            //
            //           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
            //
            Rgemv("No transpose", lastc, lastv, one, c, ldc, v, incv, zero, work, 1);
            //
            //           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
            //
            Rger(lastc, lastv, -tau, work, 1, v, incv, c, ldc);
        }
    }
    //
    //     End of Rlarf
    //
}
