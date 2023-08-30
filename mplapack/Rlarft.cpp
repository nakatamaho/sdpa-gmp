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

void Rlarft(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpf_class *v, mplapackint const ldv, mpf_class *tau, mpf_class *t, mplapackint const ldt) {
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
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    mplapackint prevlastv = 0;
    mplapackint i = 0;
    const mpf_class zero = 0.0;
    mplapackint j = 0;
    mplapackint lastv = 0;
    const mpf_class one = 1.0;
    if (Mlsame_gmp(direct, "F")) {
        prevlastv = n;
        for (i = 1; i <= k; i = i + 1) {
            prevlastv = std::max(i, prevlastv);
            if (tau[i - 1] == zero) {
                //
                //              H(i)  =  I
                //
                for (j = 1; j <= i; j = j + 1) {
                    t[(j - 1) + (i - 1) * ldt] = zero;
                }
            } else {
                //
                //              general case
                //
                if (Mlsame_gmp(storev, "C")) {
                    //                 Skip any trailing zeros.
                    for (lastv = n; lastv >= i + 1; lastv = lastv - 1) {
                        if (v[(lastv - 1) + (i - 1) * ldv] != zero) {
                            break;
                        }
                    }
                    for (j = 1; j <= i - 1; j = j + 1) {
                        t[(j - 1) + (i - 1) * ldt] = -tau[i - 1] * v[(i - 1) + (j - 1) * ldv];
                    }
                    j = std::min(lastv, prevlastv);
                    //
                    //                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
                    //
                    Rgemv("Transpose", j - i, i - 1, -tau[i - 1], &v[((i + 1) - 1)], ldv, &v[((i + 1) - 1) + (i - 1) * ldv], 1, one, &t[(i - 1) * ldt], 1);
                } else {
                    //                 Skip any trailing zeros.
                    for (lastv = n; lastv >= i + 1; lastv = lastv - 1) {
                        if (v[(i - 1) + (lastv - 1) * ldv] != zero) {
                            break;
                        }
                    }
                    for (j = 1; j <= i - 1; j = j + 1) {
                        t[(j - 1) + (i - 1) * ldt] = -tau[i - 1] * v[(j - 1) + (i - 1) * ldv];
                    }
                    j = std::min(lastv, prevlastv);
                    //
                    //                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
                    //
                    Rgemv("No transpose", i - 1, j - i, -tau[i - 1], &v[((i + 1) - 1) * ldv], ldv, &v[(i - 1) + ((i + 1) - 1) * ldv], ldv, one, &t[(i - 1) * ldt], 1);
                }
                //
                //              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
                //
                Rtrmv("Upper", "No transpose", "Non-unit", i - 1, t, ldt, &t[(i - 1) * ldt], 1);
                t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
                if (i > 1) {
                    prevlastv = std::max(prevlastv, lastv);
                } else {
                    prevlastv = lastv;
                }
            }
        }
    } else {
        prevlastv = 1;
        for (i = k; i >= 1; i = i - 1) {
            if (tau[i - 1] == zero) {
                //
                //              H(i)  =  I
                //
                for (j = i; j <= k; j = j + 1) {
                    t[(j - 1) + (i - 1) * ldt] = zero;
                }
            } else {
                //
                //              general case
                //
                if (i < k) {
                    if (Mlsame_gmp(storev, "C")) {
                        //                    Skip any leading zeros.
                        for (lastv = 1; lastv <= i - 1; lastv = lastv + 1) {
                            if (v[(lastv - 1) + (i - 1) * ldv] != zero) {
                                break;
                            }
                        }
                        for (j = i + 1; j <= k; j = j + 1) {
                            t[(j - 1) + (i - 1) * ldt] = -tau[i - 1] * v[((n - k + i) - 1) + (j - 1) * ldv];
                        }
                        j = std::max(lastv, prevlastv);
                        //
                        //                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
                        //
                        Rgemv("Transpose", n - k + i - j, k - i, -tau[i - 1], &v[(j - 1) + ((i + 1) - 1) * ldv], ldv, &v[(j - 1) + (i - 1) * ldv], 1, one, &t[((i + 1) - 1) + (i - 1) * ldt], 1);
                    } else {
                        //                    Skip any leading zeros.
                        for (lastv = 1; lastv <= i - 1; lastv = lastv + 1) {
                            if (v[(i - 1) + (lastv - 1) * ldv] != zero) {
                                break;
                            }
                        }
                        for (j = i + 1; j <= k; j = j + 1) {
                            t[(j - 1) + (i - 1) * ldt] = -tau[i - 1] * v[(j - 1) + ((n - k + i) - 1) * ldv];
                        }
                        j = std::max(lastv, prevlastv);
                        //
                        //                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
                        //
                        Rgemv("No transpose", k - i, n - k + i - j, -tau[i - 1], &v[((i + 1) - 1) + (j - 1) * ldv], ldv, &v[(i - 1) + (j - 1) * ldv], ldv, one, &t[((i + 1) - 1) + (i - 1) * ldt], 1);
                    }
                    //
                    //                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
                    //
                    Rtrmv("Lower", "No transpose", "Non-unit", k - i, &t[((i + 1) - 1) + ((i + 1) - 1) * ldt], ldt, &t[((i + 1) - 1) + (i - 1) * ldt], 1);
                    if (i > 1) {
                        prevlastv = std::min(prevlastv, lastv);
                    } else {
                        prevlastv = lastv;
                    }
                }
                t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
            }
        }
    }
    //
    //     End of Rlarft
    //
}
