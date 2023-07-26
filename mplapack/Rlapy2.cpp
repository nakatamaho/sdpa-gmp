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

#include <mpblas_dd.h>
#include <mplapack_dd.h>

dd_real Rlapy2(dd_real const x, dd_real const y) {
    dd_real return_value = 0.0;
    //
    bool x_is_nan = Risnan(x);
    bool y_is_nan = Risnan(y);
    if (x_is_nan) {
        return_value = x;
    }
    if (y_is_nan) {
        return_value = y;
    }
    //
    dd_real xabs = 0.0;
    dd_real yabs = 0.0;
    dd_real w = 0.0;
    dd_real z = 0.0;
    const dd_real zero = 0.0;
    const dd_real one = 1.0;
    if (!(x_is_nan || y_is_nan)) {
        xabs = abs(x);
        yabs = abs(y);
        w = std::max(xabs, yabs);
        z = std::min(xabs, yabs);
        if (z == zero) {
            return_value = w;
        } else {
            return_value = w * sqrt(one + pow2((z / w)));
        }
    }
    return return_value;
    //
    //     End of Rlapy2
    //
}
