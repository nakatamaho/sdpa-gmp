/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlamch.cpp,v 1.16 2010/08/07 04:48:32 nakatamaho Exp $
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
#include <stdio.h>

//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
mpf_class RlamchE_dd(void) {
    // 2^(-52-52) = 2^-104 = 4.93e-32
    return mpf_class::_eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchS_dd(void) {
    // 2^(-1022+53) = 2.0042e-292
    return mpf_class::_min_normalized;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchB_dd(void) {
    mpf_class two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchP_dd(void) {
    mpf_class base, eps, prec;

    base = RlamchB_dd();
    eps = RlamchE_dd();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchN_dd(void) {
    return (mpf_class)106.0; // 53*2
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchR_dd(void) {
    mpf_class mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchM_dd(void) { return mpf_class(-1021.0 + 53.0); }

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchU_dd(void) { return mpf_class::_min_normalized; }

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchL_dd(void) { return (mpf_class)1024.0; }

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchO_dd(void) {
    // due to bug of mpf_class, we cannot take some arithmetic for mpf_class::_max; e.g. sqrt.
    // thus we use smaller values
    mpf_class a = mpf_class::_max;
    return a;
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchZ_dd(void) {
    mpf_class mtemp = 0.0;
    return mtemp;
}

mpf_class Rlamch_dd(const char *cmach) {
    if (Mlsame_dd(cmach, "E"))
        return RlamchE_dd();
    if (Mlsame_dd(cmach, "S"))
        return RlamchS_dd();
    if (Mlsame_dd(cmach, "B"))
        return RlamchB_dd();
    if (Mlsame_dd(cmach, "P"))
        return RlamchP_dd();
    if (Mlsame_dd(cmach, "N"))
        return RlamchN_dd();
    if (Mlsame_dd(cmach, "R"))
        return RlamchR_dd();
    if (Mlsame_dd(cmach, "M"))
        return RlamchM_dd();
    if (Mlsame_dd(cmach, "U"))
        return RlamchU_dd();
    if (Mlsame_dd(cmach, "L"))
        return RlamchL_dd();
    if (Mlsame_dd(cmach, "O"))
        return RlamchO_dd();

    Mxerbla_dd("Rlamch", 1);
    return RlamchZ_dd();
}

mpf_class Rlamc3(mpf_class a, mpf_class b) { return a + b; }
