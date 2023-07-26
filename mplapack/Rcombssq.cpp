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

void Rcombssq(mpf_class *v1, mpf_class *v2) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //     November 2018
    //
    //     .. Array Arguments ..
    //     ..
    //
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Executable Statements ..
    //
    const mpf_class zero = 0.0;
    if (v1[1 - 1] >= v2[1 - 1]) {
        if (v1[1 - 1] != zero) {
            v1[2 - 1] += pow2((v2[1 - 1] / v1[1 - 1])) * v2[2 - 1];
        } else {
            v1[2 - 1] += v2[2 - 1];
        }
    } else {
        v1[2 - 1] = v2[2 - 1] + pow2((v1[1 - 1] / v2[1 - 1])) * v1[2 - 1];
        v1[1 - 1] = v2[1 - 1];
    }
    //
    //     End of Rcombssq
    //
}
