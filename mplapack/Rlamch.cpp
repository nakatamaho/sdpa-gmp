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
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>
#include <stdio.h>

//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
mpf_class RlamchE_gmp(void) {
    mpf_class eps;
    mpf_class one;
    unsigned long exp2;
    one = 1.0;
    exp2 = mpf_get_prec(one.get_mpf_t());
    mpf_div_2exp(eps.get_mpf_t(), one.get_mpf_t(), exp2);
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchS_gmp(void) {
    mpf_class sfmin;
    mpf_class one = 1.0;
    unsigned long exp2;
#if defined _WIN32 //exponent investigated by MPFR
    exp2 = 0x3FFFFFFF;
#else
    exp2 = 0x3FFFFFFFFFFFFFFFL;
#endif
    mpf_div_2exp(sfmin.get_mpf_t(), one.get_mpf_t(), exp2);
    return sfmin;

    /* following code fragment is to test safe minimum
        mpf_class largenum;
        for(int p = 60; p>=0; p--) {
          for (int a=16; a<= 5120; a=a+128) {
            sfmin = 0.0;
            mpf_set_default_prec(a);
            exp2 = (1UL << (mp_bits_per_limb-p)) - 1;
            mpf_div_2exp(sfmin.get_mpf_t(), one.get_mpf_t(), exp2);
            largenum  = 1.0 / sfmin;
            gmp_printf("%d, a:%5d, p:%5d sfmin: %16.20Fe largenum: %16.20Fe\n", mp_bits_per_limb, a, p, sfmin.get_mpf_t(), largenum.get_mpf_t());
            if (sfmin < 1.0 ) { printf("sfmin yes\n"); }
            if (largenum > 1.0 ) { printf("largenum yes\n"); }
          }
        }
    */
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchB_gmp(void) {
    mpf_class two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchP_gmp(void) {
    mpf_class base, eps, prec;

    base = RlamchB_gmp();
    eps = RlamchE_gmp();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchN_gmp(void) {
    unsigned long int tmp;
    mpf_class mtmp;
    mpf_class mtmp2;
    tmp = mpf_get_prec(mtmp.get_mpf_t());
    mtmp2 = tmp;
    return mtmp2;

    /* following is fragment of code to test digits in mantissa
       for (int a=8; a<= 5120; a=a+128) {
         mpf_set_default_prec(a);
         mpf_class mmtmp;
         tmp = mpf_get_prec(mmtmp.get_mpf_t());
         printf("number of digits in mantissa %d\n", (int)tmp );
       }
    */
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchR_gmp(void) {
    // always rounding in addition on GMP.
    mpf_class mtmp;

    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchM_gmp(void) {
    unsigned long exp2;
    mpf_class tmp;
    mpf_class uflowmin, one = 1.0;
#if defined _WIN32 //exponent investigated by MPFR
    exp2 = 0x3FFFFFFF;
#else
    exp2 = 0x3FFFFFFFFFFFFFFFL; //exponent investigated by MPFR
#endif
    tmp = exp2;
    return -tmp;

    /*
      following code fragment is to test minimum
      exponent before (gradual) underflow...but we just got Bus error
        mpf_class mtmp1, mtmp2;
        for(int p = 11; p>=0; p--) {
          for (int a=51200; a<= 102400; a=a+128) {
            mpf_set_default_prec(a);
            printf("p %d a %d \n", p, a);
            uflowmin=0.0;
            exp2 = (1UL << (mp_bits_per_limb-p)) - 1;
            mpf_div_2exp(uflowmin.get_mpf_t(), one.get_mpf_t(), exp2);
            mtmp1 = uflowmin/2.0;
            gmp_printf("p %d, uflowmin: %16.20Fe uflowmin/2 %16.20Fe\n", p, uflowmin.get_mpf_t(), mtmp1.get_mpf_t());
            mtmp2 = mtmp1 * 2.0;
            gmp_printf("mtmp2: %16.20Fe %lu\n", mtmp2.get_mpf_t(),exp2);
            if (uflowmin != mtmp2 ) { printf("underflow\n"); exit(1); }
          }
        }
    */
}

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchU_gmp(void) {
    mpf_class underflowmin;
    mpf_class one = 1.0;
    unsigned long exp2;
#if defined _WIN32
    exp2 = 0x3FFFFFFF;
#else
    exp2 = 0x3FFFFFFFFFFFFFFFL;
#endif
    mpf_div_2exp(underflowmin.get_mpf_t(), one.get_mpf_t(), exp2);
    return underflowmin;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchL_gmp(void) {
    mpf_class maxexp;
    unsigned long exp2;
#if defined _WIN32
    exp2 = 0x3FFFFFFF;
#else
    exp2 = 0x3FFFFFFFFFFFFFFFL;
#endif
    return maxexp;
}

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchO_gmp(void) {
    mpf_class overflowmax;
    mpf_class one = 1.0;
    unsigned long exp2;
#if defined _WIN32
    exp2 = 0x3FFFFFFF;
#else
    exp2 = 0x3FFFFFFFFFFFFFFFL;
#endif
    mpf_mul_2exp(overflowmax.get_mpf_t(), one.get_mpf_t(), exp2);
    return overflowmax;
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
mpf_class RlamchZ_gmp(void) {
    mpf_class mtemp = 0.0;
    return mtemp;
}

mpf_class Rlamch_gmp(const char *cmach) {
    if (Mlsame_gmp(cmach, "E"))
        return RlamchE_gmp();
    if (Mlsame_gmp(cmach, "S"))
        return RlamchS_gmp();
    if (Mlsame_gmp(cmach, "B"))
        return RlamchB_gmp();
    if (Mlsame_gmp(cmach, "P"))
        return RlamchP_gmp();
    if (Mlsame_gmp(cmach, "N"))
        return RlamchN_gmp();
    if (Mlsame_gmp(cmach, "R"))
        return RlamchR_gmp();
    if (Mlsame_gmp(cmach, "M"))
        return RlamchM_gmp();
    if (Mlsame_gmp(cmach, "U"))
        return RlamchU_gmp();
    if (Mlsame_gmp(cmach, "L"))
        return RlamchL_gmp();
    if (Mlsame_gmp(cmach, "O"))
        return RlamchO_gmp();

    Mxerbla_gmp("Rlamch", 1);
    return RlamchZ_gmp();
}

mpf_class Rlamc3(mpf_class a, mpf_class b) { return a + b; }
