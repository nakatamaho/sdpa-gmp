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

mpf_class fast_approx_log(const mpf_class& value) {
    // Get the mantissa (a) and exponent (b) such that value = a * 2^b
    mp_exp_t exponent;
    double mantissa = mpf_get_d_2exp(&exponent, value.get_mpf_t());

    // Calculate log(a) + b * log(2)
    double log_mantissa = std::log(mantissa);
    double log_2 = std::log(2.0);
    return log_mantissa + exponent * log_2;
}

void Rlartg(mpf_class const f, mpf_class const g, mpf_class &cs, mpf_class &sn, mpf_class &r) {
    mpf_class safmin = 0.0;
    mpf_class eps = 0.0;
    const mpf_class two = 2.0;
    mpf_class safmn2 = 0.0;
    const mpf_class one = 1.0;
    mpf_class safmx2 = 0.0;
    const mpf_class zero = 0.0;
    mpf_class f1 = 0.0;
    mpf_class g1 = 0.0;
    mpf_class scale = 0.0;
    mplapackint count = 0;
    mplapackint i = 0;
    //
    safmin = Rlamch_gmp("S");
    eps = Rlamch_gmp("E");

    //mpf_class safmn2_org = gmpxx::pow(Rlamch_gmp("B"), castINTEGER_gmp(gmpxx::log(safmin / eps) / gmpxx::log(Rlamch_gmp("B")) / two));

    mpf_class tmp = -fast_approx_log(safmin / eps) / fast_approx_log(Rlamch_gmp("B")) / two;
    safmn2 = 1;
    //safmn2.div_2exp(castINTEGER_gmp(tmp)); //div_2exp accepts unsigned long
    mpf_div_2exp(safmn2.get_mpf_t(), safmn2.get_mpf_t(), castINTEGER_gmp(tmp));

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
