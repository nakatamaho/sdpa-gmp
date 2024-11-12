/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_gmp.h,v 1.10 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_GMP_H_
#define _MPBLAS_GMP_H_

#define ___MPLAPACK_GMP_DEFAULT_PRECISION___ 512

#include "gmpxx.h"
#include "mplapack_config.h"
#include "mplapack_utils_gmp.h"

mpf_class Rdot(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Rcopy(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Raxpy(mplapackint const n, mpf_class const da, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Rscal(mplapackint const n, mpf_class const da, mpf_class *dx, mplapackint const incx);
bool Mlsame_gmp(const char *a, const char *b);
void Mxerbla_gmp(const char *srname, int info);
void Rswap(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
mpf_class Rnrm2(mplapackint const n, mpf_class *x, mplapackint const incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rsymv(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rsyr2(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *y, mplapackint const incy, mpf_class *a, mplapackint const lda);
void Rger(mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *y, mplapackint const incy, mpf_class *a, mplapackint const lda);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc);

#endif
