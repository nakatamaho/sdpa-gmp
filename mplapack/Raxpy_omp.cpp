/*
 * Copyright (c) 2010-2025
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Raxpy.cpp,v 1.11 2010/08/07 05:50:09 nakatamaho Exp $
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

/*
Based on http://www.netlib.org/blas/daxpy.f
*/

#include <mpblas_gmp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void Raxpy_omp(mplapackint n, mpf_class da, mpf_class *dx, mplapackint incx, mpf_class *dy, mplapackint incy) {
    mplapackint ix = 0;
    mplapackint iy = 0;

    if (n <= 0)
        return;
    if (da == 0.0)
        return;

    if (incx < 0)
        ix = (-n + 1) * incx;
    if (incy < 0)
        iy = (-n + 1) * incy;

    if (incx == 1 && incy == 1) {
#ifdef _OPENMP
//#pragma omp parallel
#endif
        {
#ifdef _OPENMP
//#pragma omp for
#endif
            for (mplapackint i = 0; i < n; i++) {
                mpf_class temp = dx[i];
                temp *= da;
                dy[i] += temp;
            }
        }
    } else {
        for (mplapackint i = 0; i < n; i++) {
            mpf_class temp = dx[ix];
            temp *= da;
            dy[iy] += temp;
            ix += incx;
            iy += incy;
        }
    }
}
