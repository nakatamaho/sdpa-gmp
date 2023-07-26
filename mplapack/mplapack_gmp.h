/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_dd.h,v 1.31 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_GMP_H_
#define _MPLAPACK_GMP_H_

#include "mplapack_config.h"
#include "gmpxx.h"

/* this is a subset of mplapack only for SDPA-GMP */

mpf_class Rlamch_dd(const char *cmach);
void Rsyev(const char *jobz, const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *w, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rsteqr(const char *compz, mplapackint const n, mpf_class *d, mpf_class *e, mpf_class *z, mplapackint const ldz, mpf_class *work, mplapackint &info);
void Rpotrf(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info);
void Rpotf2(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info);
void Rpotrf2(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info);
void Rlascl(const char *type, mplapackint const kl, mplapackint const ku, mpf_class const cfrom, mpf_class const cto, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info);
void Rlascl2(mplapackint const m, mplapackint const n, mpf_class *d, mpf_class *x, mplapackint const ldx);
void Rsytrd(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *d, mpf_class *e, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rsytrd_2stage(const char *vect, const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *d, mpf_class *e, mpf_class *tau, mpf_class *hous2, mplapackint const lhous2, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rsytrd_sb2st(const char *stage1, const char *vect, const char *uplo, mplapackint const n, mplapackint const kd, mpf_class *ab, mplapackint const ldab, mpf_class *d, mpf_class *e, mpf_class *hous, mplapackint const lhous, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rsytrd_sy2sb(const char *uplo, mplapackint const n, mplapackint const kd, mpf_class *a, mplapackint const lda, mpf_class *ab, mplapackint const ldab, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rsytd2(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *d, mpf_class *e, mpf_class *tau, mplapackint &info);
void Rlae2(mpf_class const a, mpf_class const b, mpf_class const c, mpf_class &rt1, mpf_class &rt2);
void Rlasrt(const char *id, mplapackint const n, mpf_class *d, mplapackint &info);
void Rorgql(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rorgqr(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info);
void Rlarfg(mplapackint const n, mpf_class &alpha, mpf_class *x, mplapackint const incx, mpf_class &tau);
void Rlarfgp(mplapackint const n, mpf_class &alpha, mpf_class *x, mplapackint const incx, mpf_class &tau);
void Rlassq(mplapackint const n, mpf_class *x, mplapackint const incx, mpf_class &scale, mpf_class &sumsq);
void Rorg2l(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint &info);
void Rlarft(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpf_class *v, mplapackint const ldv, mpf_class *tau, mpf_class *t, mplapackint const ldt);
void Rlarfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *v, mplapackint const ldv, mpf_class *t, mplapackint const ldt, mpf_class *c, mplapackint const ldc, mpf_class *work, mplapackint const ldwork);
void Rlarfb_gett(const char *ident, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *t, mplapackint const ldt, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class *work, mplapackint const ldwork);
void Rorg2r(mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint &info);
void Rlarf(const char *side, mplapackint const m, mplapackint const n, mpf_class *v, mplapackint const incv, mpf_class const tau, mpf_class *c, mplapackint const ldc, mpf_class *work);
void Rlarfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *v, mplapackint const ldv, mpf_class *t, mplapackint const ldt, mpf_class *c, mplapackint const ldc, mpf_class *work, mplapackint const ldwork);
void Rlarfb_gett(const char *ident, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class *t, mplapackint const ldt, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class *work, mplapackint const ldwork);
void Rlarfg(mplapackint const n, mpf_class &alpha, mpf_class *x, mplapackint const incx, mpf_class &tau);
void Rlarfgp(mplapackint const n, mpf_class &alpha, mpf_class *x, mplapackint const incx, mpf_class &tau);
void Rlarft(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpf_class *v, mplapackint const ldv, mpf_class *tau, mpf_class *t, mplapackint const ldt);
void Rlarfx(const char *side, mplapackint const m, mplapackint const n, mpf_class *v, mpf_class const tau, mpf_class *c, mplapackint const ldc, mpf_class *work);
void Rlarfy(const char *uplo, mplapackint const n, mpf_class *v, mplapackint const incv, mpf_class const tau, mpf_class *c, mplapackint const ldc, mpf_class *work);
void Rpotf2(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint &info);
void Rlaset(const char *uplo, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class const beta, mpf_class *a, mplapackint const lda);
void Rlaev2(mpf_class const a, mpf_class const b, mpf_class const c, mpf_class &rt1, mpf_class &rt2, mpf_class &cs1, mpf_class &sn1);
void Rlasr(const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpf_class *c, mpf_class *s, mpf_class *a, mplapackint const lda);
void Rlasrt(const char *id, mplapackint const n, mpf_class *d, mplapackint &info);
void Rlartg(mpf_class const f, mpf_class const g, mpf_class &cs, mpf_class &sn, mpf_class &r);
void Rlartgp(mpf_class const f, mpf_class const g, mpf_class &cs, mpf_class &sn, mpf_class &r);
void Rlartgs(mpf_class const x, mpf_class const y, mpf_class const sigma, mpf_class &cs, mpf_class &sn);
void Rlatrd(const char *uplo, mplapackint const n, mplapackint const nb, mpf_class *a, mplapackint const lda, mpf_class *e, mpf_class *tau, mpf_class *w, mplapackint const ldw);
void Rsterf(mplapackint const n, mpf_class *d, mpf_class *e, mplapackint &info);
void Rorgtr(const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *tau, mpf_class *work, mplapackint const lwork, mplapackint &info);
mplapackint iMparmq_dd(mplapackint const ispec, const char *name, const char *opts, mplapackint const n, mplapackint const ilo, mplapackint const ihi, mplapackint const lwork);
mplapackint iMieeeck_dd(mplapackint const &ispec, mpf_class const &zero, mpf_class const &one);
bool Rlaisnan(mpf_class const din1, mpf_class const din2);
bool Risnan(mpf_class const din);
mplapackint iMlaenv_dd(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mpf_class Rlanst(const char *norm, mplapackint const n, mpf_class *d, mpf_class *e);
mpf_class Rlapy2(mpf_class const x, mpf_class const y);
mplapackint iMladlr(mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda);
mplapackint iMladlc(mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda);
mpf_class Rlansy(const char *norm, const char *uplo, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *work);
void Rcombssq(mpf_class *v1, mpf_class *v2);

#endif
