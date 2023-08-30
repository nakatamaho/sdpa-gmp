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

void Rlasr(const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpf_class *c, mpf_class *s, mpf_class *a, mplapackint const lda) {
    //
    //     Test the input parameters
    //
    mplapackint info = 0;
    if (!(Mlsame_gmp(side, "L") || Mlsame_gmp(side, "R"))) {
        info = 1;
    } else if (!(Mlsame_gmp(pivot, "V") || Mlsame_gmp(pivot, "T") || Mlsame_gmp(pivot, "B"))) {
        info = 2;
    } else if (!(Mlsame_gmp(direct, "F") || Mlsame_gmp(direct, "B"))) {
        info = 3;
    } else if (m < 0) {
        info = 4;
    } else if (n < 0) {
        info = 5;
    } else if (lda < std::max((mplapackint)1, m)) {
        info = 9;
    }
    if (info != 0) {
        Mxerbla_gmp("Rlasr", info);
        return;
    }
    //
    //     Quick return if possible
    //
    if ((m == 0) || (n == 0)) {
        return;
    }
    mplapackint j = 0;
    mpf_class ctemp = 0.0;
    mpf_class stemp = 0.0;
    const mpf_class one = 1.0;
    const mpf_class zero = 0.0;
    mplapackint i = 0;
    mpf_class temp = 0.0;
    if (Mlsame_gmp(side, "L")) {
        //
        //        Form  P * A
        //
        if (Mlsame_gmp(pivot, "V")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 1; j <= m - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[((j + 1) - 1) + (i - 1) * lda];
                            a[((j + 1) - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * temp + ctemp * a[(j - 1) + (i - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = m - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[((j + 1) - 1) + (i - 1) * lda];
                            a[((j + 1) - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * temp + ctemp * a[(j - 1) + (i - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame_gmp(pivot, "T")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 2; j <= m; j = j + 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(i - 1) * lda];
                            a[(i - 1) * lda] = stemp * temp + ctemp * a[(i - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = m; j >= 2; j = j - 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(i - 1) * lda];
                            a[(i - 1) * lda] = stemp * temp + ctemp * a[(i - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame_gmp(pivot, "B")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 1; j <= m - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * a[(m - 1) + (i - 1) * lda] + ctemp * temp;
                            a[(m - 1) + (i - 1) * lda] = ctemp * a[(m - 1) + (i - 1) * lda] - stemp * temp;
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = m - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * a[(m - 1) + (i - 1) * lda] + ctemp * temp;
                            a[(m - 1) + (i - 1) * lda] = ctemp * a[(m - 1) + (i - 1) * lda] - stemp * temp;
                        }
                    }
                }
            }
        }
    } else if (Mlsame_gmp(side, "R")) {
        //
        //        Form A * P**T
        //
        if (Mlsame_gmp(pivot, "V")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + ((j + 1) - 1) * lda];
                            a[(i - 1) + ((j + 1) - 1) * lda] = ctemp * temp - stemp * a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * temp + ctemp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = n - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + ((j + 1) - 1) * lda];
                            a[(i - 1) + ((j + 1) - 1) * lda] = ctemp * temp - stemp * a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * temp + ctemp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame_gmp(pivot, "T")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 2; j <= n; j = j + 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = ctemp * temp - stemp * a[(i - 1)];
                            a[(i - 1)] = stemp * temp + ctemp * a[(i - 1)];
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = n; j >= 2; j = j - 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = ctemp * temp - stemp * a[(i - 1)];
                            a[(i - 1)] = stemp * temp + ctemp * a[(i - 1)];
                        }
                    }
                }
            }
        } else if (Mlsame_gmp(pivot, "B")) {
            if (Mlsame_gmp(direct, "F")) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * a[(i - 1) + (n - 1) * lda] + ctemp * temp;
                            a[(i - 1) + (n - 1) * lda] = ctemp * a[(i - 1) + (n - 1) * lda] - stemp * temp;
                        }
                    }
                }
            } else if (Mlsame_gmp(direct, "B")) {
                for (j = n - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * a[(i - 1) + (n - 1) * lda] + ctemp * temp;
                            a[(i - 1) + (n - 1) * lda] = ctemp * a[(i - 1) + (n - 1) * lda] - stemp * temp;
                        }
                    }
                }
            }
        }
    }
    //
    //     End of Rlasr
    //
}
