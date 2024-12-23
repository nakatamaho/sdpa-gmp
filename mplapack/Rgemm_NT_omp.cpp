/*
 * Copyright (c) 2008-2024
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm_NT.cpp,v 1.1 2010/12/28 06:13:53 nakatamaho Exp $
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
#ifdef _OPENMP
#include <omp.h>
#endif

void Rgemm_NT_omp(mplapackint m, mplapackint n, mplapackint k, mpf_class alpha, mpf_class *A, mplapackint lda, mpf_class *B, mplapackint ldb, mpf_class beta, mpf_class *C, mplapackint ldc) {
    // Form C := alpha*A*B' + beta*C.
    mplapackint i, j, l;
    mpf_class temp, templ;

    // Scale C by beta: C = alpha * A * B' + beta * C
    // Handle special cases for beta to optimize performance
    if (beta == 0.0) {
        // If beta is 0, set C to zero
//#pragma omp parallel for collapse(2) schedule(static) private(i, j)
        for (j = 0; j < n; j++) {
            for (i = 0; i < m; i++) {
                C[i + j * ldc] = 0.0;
            }
        }
    } else if (beta != 1.0) {
        // If beta is not 1, scale C by beta
//#pragma omp parallel for collapse(2) schedule(static) private(i, j)
        for (j = 0; j < n; j++) {
            for (i = 0; i < m; i++) {
                C[i + j * ldc] *= beta;
            }
        }
    }
    // If beta is 1, no scaling is needed

// Compute alpha * A * B' and add to C: C += alpha * A * B'
//#pragma omp parallel for private(j, l, i, temp, templ) schedule(static)
    for (j = 0; j < n; j++) {
        for (l = 0; l < k; l++) {
            temp = alpha;
            temp *= B[j + l * ldb];
            for (i = 0; i < m; i++) {
                templ = temp;
                templ *= A[i + l * lda];
                C[i + j * ldc] += templ;
            }
        }
    }
    return;
}
