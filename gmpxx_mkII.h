/*
 * Copyright (c) 2024
 *      Nakata, Maho
 *      All rights reserved.
 *
 *
 * The gmpxx_mkII.h is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at your
 * option) any later version.
 *
 * The gmpxx_mkII.h is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the gmpxx_mkII.h; see the file LICENSE.  If not, see
 * http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef ___GMPXX_MKII_H___
#define ___GMPXX_MKII_H___

#if __cplusplus < 201703L
#error "This class only runs on C++ 17 and later"
#endif

#include <gmp.h>
#include <limits>
#include <iostream>
#include <utility>
#include <cassert>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <tuple>
#include <iomanip>

#define ___MPF_CLASS_EXPLICIT___ explicit

#if defined ___GMPXX_STRICT_COMPATIBILITY___
#define ___GMPXX_DONT_USE_NAMESPACE___
#define ___GMPXX_UDL_CHAR___
#endif

#define INT_COND(T, X) typename std::enable_if<std::is_integral<T>::value, X>::type
#define FLOAT_COND(T, X) typename std::enable_if<std::is_float<T>::value, X>::type
#define UNSIGNED_INT_COND(T, X) typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, X>::type
#define SIGNED_INT_COND(T, X) typename std::enable_if<std::is_integral<T>::value && std::is_signed<T>::value, X>::type
#define NON_INT_COND(T, X) typename std::enable_if<std::is_arithmetic<T>::value && !std::is_integral<T>::value, X>::type
#define NON_GMP_COND(T, X) typename std::enable_if<!std::is_same<T, mpf_class>::value && !std::is_same<T, mpq_class>::value && !std::is_same<T, mpz_class>::value, X>::type

#if !defined ___GMPXX_DONT_USE_NAMESPACE___
namespace gmpxx {
#endif

class mpz_class;
class mpq_class;
class mpf_class;

struct gmpxx_defaults {
    static void set_default_prec(int prec) { mpf_set_default_prec(prec); }
    static mp_bitcnt_t get_default_prec() { return mpf_get_default_prec(); }
    inline static int base = 10;
};
class mpf_class_initializer {
  public:
    mpf_class_initializer() {
        int prec = 512;

        const char *prec_env = std::getenv("GMPXX_MKII_DEFAULT_PREC");

        if (prec_env) {
            if (is_positive_integer(prec_env)) {
                int prec_val = std::stoi(prec_env);
                if (prec_val > 0) {
                    prec = prec_val;
                }
            } else {
                std::cerr << "Error: Invalid GMPXX_MKII_DEFAULT_PREC value: must be a positive integer." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }

        mpf_set_default_prec(prec);
        gmpxx_defaults::base = 10;
    }

  private:
    bool is_positive_integer(const std::string &s) { return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit); }
};
inline mpf_class_initializer global_mpf_class_initializer;
class mpf_class_initializer_singleton {
  public:
    static mpf_class_initializer_singleton &instance() {
        static mpf_class_initializer_singleton instance;
        return instance;
    }
    static void initialize() { instance(); }

  private:
    mpf_class_initializer_singleton() { (void)global_mpf_class_initializer; }
    mpf_class_initializer_singleton(const mpf_class_initializer_singleton &) = delete;
    mpf_class_initializer_singleton &operator=(const mpf_class_initializer_singleton &) = delete;
};
namespace {
const auto &initializer = mpf_class_initializer_singleton::instance();
}
#define ___GMPXX_MKII_INITIALIZER___ mpf_class_initializer_singleton::instance()

template <typename T = void> struct caches {
    static mpf_class pi_cached;
    static mpf_class e_cached;
    static mpf_class log_cached;
    static mpf_class log2_cached;
};
template <typename T> mpf_class caches<T>::pi_cached;
template <typename T> mpf_class caches<T>::e_cached;
template <typename T> mpf_class caches<T>::log_cached;
template <typename T> mpf_class caches<T>::log2_cached;

class mpz_class {
  public:
    ////////////////////////////////////////////////////////////////////////////////////////
    // 12.2 C++ Interface Integers
    // cf. https://gmplib.org/manual/C_002b_002b-Interface-Integers
    ////////////////////////////////////////////////////////////////////////////////////////
    // constructors and destructors
    mpz_class() { mpz_init(value); }
    // The rule of 0/3/5
    // The rule 1 of 5 copy constructor
    mpz_class(const mpz_class &op) {
        mpz_init(value);
        mpz_set(value, op.value);
    }
    // The rule 2 of 5 copy assignment operator
    mpz_class &operator=(const mpz_class &op) noexcept {
        if (this != &op) {
            mpz_set(value, op.value);
        }
        return *this;
    }
    // The rule 3 of 5 default deconstructor
    ~mpz_class() { mpz_clear(value); }
    // The rule 4 of 5 move constructor
    mpz_class(mpz_class &&op) noexcept {
        mpz_init(value);
        mpz_swap(value, op.value);
    }
    // The rule 5 of 5 move assignment operator
    mpz_class &operator=(mpz_class &&op) noexcept {
        if (this != &op) {
            mpz_swap(value, op.value);
        }
        return *this;
    }
    // constructors
    explicit mpz_class(const mpz_t z) {
        mpz_init(value);
        mpz_set(value, z);
    }
    mpz_class(const mpq_t op) {
        mpz_init(value);
        mpz_set_q(value, op);
    }
    mpz_class(const mpf_t op) {
        mpz_init(value);
        mpz_set_f(value, op);
    }
    mpz_class(const char *str, int base = 0) {
        mpz_init(value);
        if (mpz_set_str(value, str, base) != 0) {
            throw std::invalid_argument("");
        }
    }
    mpz_class(const std::string &str, int base = 0) {
        mpz_init(value);
        if (mpz_set_str(value, str.c_str(), base) != 0) {
            throw std::invalid_argument("");
        }
    }
    mpz_class(unsigned long int op) { mpz_init_set_ui(value, op); }
    mpz_class(signed long int op) { mpz_init_set_si(value, op); }
    mpz_class(double op) { mpz_init_set_d(value, op); }
    mpz_class(unsigned int op) { mpz_init_set_ui(value, static_cast<unsigned long int>(op)); }
    mpz_class(signed int op) { mpz_init_set_si(value, static_cast<signed long int>(op)); }
    // assignments from other objects
    mpz_class &operator=(double d) noexcept {
        mpz_set_d(value, d);
        return *this;
    }
    mpz_class &operator=(const char *str) {
        if (mpz_set_str(value, str, 0) != 0) {
            throw std::invalid_argument("");
        }
        return *this;
    }
    mpz_class &operator=(const std::string &str) {
        if (mpz_set_str(value, str.c_str(), 0) != 0) {
            throw std::invalid_argument("");
        }
        return *this;
    }
    mpz_class &operator=(const signed long int op);
    mpz_class &operator=(const unsigned long int op);
    mpz_class &operator=(const signed int op);
    mpz_class &operator=(const unsigned int op);
    mpz_class &operator=(const signed char op);
    mpz_class &operator=(const unsigned char op);
    mpz_class &operator=(const char op);
    // operators
    mpz_class operator~() const {
        mpz_class result;
        mpz_com(result.value, value);
        return result;
    }
    mpz_class &operator++() {
        mpz_add_ui(value, value, 1);
        return *this;
    }
    mpz_class operator++(int) {
        mpz_class original = *this;
        ++(*this);
        return original;
    }
    mpz_class &operator--() {
        mpz_sub_ui(value, value, 1);
        return *this;
    }
    mpz_class operator--(int) {
        mpz_class original = *this;
        --(*this);
        return original;
    }
    template <typename T> INT_COND(T, mpz_class &) operator<<=(T n) {
        mpz_mul_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> INT_COND(T, mpz_class &) operator>>=(T n) {
        mpz_tdiv_q_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> friend INT_COND(T, mpz_class) operator<<(const mpz_class &op1, T op2) {
        mpz_class result(op1);
        mpz_mul_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }
    template <typename T> friend INT_COND(T, mpz_class) operator>>(const mpz_class &op1, T op2) {
        mpz_class result(op1);
        mpz_fdiv_q_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }

    // mpz_class arithmetic operators
    inline friend mpz_class &operator+=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator-=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator*=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator/=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator%=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator&=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator|=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class &operator^=(mpz_class &lhs, const mpz_class &rhs);
    inline friend mpz_class operator+(const mpz_class &op);
    inline friend mpz_class operator-(const mpz_class &op);
    inline friend mpz_class operator+(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator-(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator*(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator/(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator%(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator&(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator|(const mpz_class &op1, const mpz_class &op2);
    inline friend mpz_class operator^(const mpz_class &op1, const mpz_class &op2);

    // mpz_class comparison operators
    inline friend bool operator==(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) == 0; }
    inline friend bool operator!=(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) != 0; }
    inline friend bool operator<(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) < 0; }
    inline friend bool operator>(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) > 0; }
    inline friend bool operator<=(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) <= 0; }
    inline friend bool operator>=(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.value, op2.value) >= 0; }

    inline friend bool operator==(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) == 0; }
    inline friend bool operator!=(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) != 0; }
    inline friend bool operator<(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) < 0; }
    inline friend bool operator>(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) > 0; }
    inline friend bool operator<=(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) <= 0; }
    inline friend bool operator>=(const mpz_class &op1, const double &op2) { return mpz_cmp_d(op1.value, op2) >= 0; }

    inline friend bool operator==(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) == 0; }
    inline friend bool operator!=(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) != 0; }
    inline friend bool operator<(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) > 0; }
    inline friend bool operator>(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) < 0; }
    inline friend bool operator<=(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) >= 0; }
    inline friend bool operator>=(double &op1, const mpz_class &op2) { return mpz_cmp_d(op2.value, op1) <= 0; }

    // mpz_class comparison operators (template version)
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator==(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) == 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator!=(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) != 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) < 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) > 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<=(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) <= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>=(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) >= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator==(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) == 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator!=(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) != 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) > 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) < 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<=(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) >= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>=(T op1, const mpz_class &op2) { return mpz_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) <= 0; }

    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator==(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) == 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator!=(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) != 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) < 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) > 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<=(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) <= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>=(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.value, static_cast<signed long int>(op2)) >= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator==(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) == 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator!=(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) != 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) > 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) < 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<=(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) >= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>=(T op1, const mpz_class &op2) { return mpz_cmp_si(op2.value, static_cast<signed long int>(op1)) <= 0; }

    template <typename T> inline friend NON_INT_COND(T, bool) operator==(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) == 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator!=(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) != 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator<(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) < 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator>(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) > 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator<=(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) <= 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator>=(const mpz_class &op1, T op2) { return mpz_cmp(op1.value, mpz_class(op2).get_mpz_t()) >= 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator==(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) == 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator!=(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) != 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator<(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) > 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator>(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) < 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator<=(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) >= 0; }
    template <typename T> inline friend NON_INT_COND(T, bool) operator>=(T op1, const mpz_class &op2) { return mpz_cmp(op2.value, mpz_class(op1).get_mpz_t()) <= 0; }

    // mpz_class arithmetic and logical operators (template version)
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2);

    // mpz_class abs (mpz_class op)
    inline friend mpz_class abs(const mpz_class &op);

    // bool mpz_class::fits_sint_p (void)
    // bool mpz_class::fits_slong_p (void)
    // bool mpz_class::fits_sshort_p (void)
    // bool mpz_class::fits_uint_p (void)
    // bool mpz_class::fits_ulong_p (void)
    // bool mpz_class::fits_ushort_p (void)
    inline bool fits_sint_p() const { return mpz_fits_sint_p(value) != 0; }
    inline bool fits_slong_p() const { return mpz_fits_slong_p(value) != 0; }
    inline bool fits_sshort_p() const { return mpz_fits_sshort_p(value) != 0; }
    inline bool fits_uint_p() const { return mpz_fits_uint_p(value) != 0; }
    inline bool fits_ulong_p() const { return mpz_fits_ulong_p(value) != 0; }
    inline bool fits_ushort_p() const { return mpz_fits_ushort_p(value) != 0; }

    // double mpz_class::get_d (void)
    // long mpz_class::get_si (void)
    // unsigned long mpz_class::get_ui (void)
    inline double get_d() const { return mpz_get_d(value); }
    inline signed long int get_si() const { return mpz_get_si(value); }
    inline unsigned long int get_ui() const { return mpz_get_ui(value); }
    // string mpz_class::get_str (int base = 10)
    inline std::string get_str(int base = 10) const {
        char *temp = mpz_get_str(nullptr, base, value);
        std::string result(temp);
        void (*freefunc)(void *, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(temp, std::strlen(temp) + 1);
        return result;
    }
    // int mpz_class::set_str (const char *str, int base)
    // int mpz_class::set_str (const string& str, int base)
    int set_str(const char *str, int base) { return mpz_set_str(value, str, base); }
    int set_str(const std::string &str, int base) { return mpz_set_str(value, str.c_str(), base); }

    // int sgn (mpz_class op)
    // mpz_class sqrt (mpz_class op)
    // mpz_class gcd (mpz_class op1, mpz_class op2)
    // mpz_class lcm (mpz_class op1, mpz_class op2)
    friend int sgn(const mpz_class &op);
    friend mpz_class sqrt(const mpz_class &op);
    friend mpz_class gcd(const mpz_class &op1, const mpz_class &op2);
    friend mpz_class lcm(const mpz_class &op1, const mpz_class &op2);

    // mpz_class mpz_class::factorial (type op)
    // mpz_class factorial (mpz_class op)
    // mpz_class mpz_class::fibonacci (type op)
    // mpz_class fibonacci (mpz_class op)
    // mpz_class mpz_class::primorial (type op)
    // mpz_class primorial (mpz_class op)
    static mpz_class factorial(const mpz_class &n) {
        if (n < 0) {
            throw std::domain_error("factorial(negative)");
        }
        double log2_n = mpz_sizeinbase(n.get_mpz_t(), 2);
        if (log2_n > 300) {
            throw std::bad_alloc();
        }
        mpz_class result;
        try {
            mpz_fac_ui(result.get_mpz_t(), n.get_ui());
        } catch (const std::bad_alloc &) {
            throw;
        }
        return result;
    }
    static mpz_class primorial(const mpz_class &op) {
        if (op < 0) {
            throw std::domain_error("primorial(negative)");
        }
        double log2_n = mpz_sizeinbase(op.get_mpz_t(), 2);
        if (log2_n > 300) {
            throw std::bad_alloc();
        }
        mpz_class result;
        try {
            mpz_primorial_ui(result.get_mpz_t(), op.get_ui());
        } catch (const std::bad_alloc &) {
            throw;
        }
        return result;
    }
    friend mpz_class primorial(const mpz_class &op);
    static mpz_class fibonacci(const mpz_class &op) {
        double log2_op = mpz_sizeinbase(op.get_mpz_t(), 2);
        if (log2_op > 300) {
            throw std::bad_alloc();
        }
        mpz_class adjusted_op = op;
        bool isNegative = op < 0;
        if (isNegative) {
            adjusted_op = -op;
        }
        unsigned long int n = adjusted_op.get_ui();
        mpz_class result;
        mpz_fib_ui(result.get_mpz_t(), n);
        if (isNegative) {
            if ((op + 1) % 2 != 0) {
                result = -result;
            }
        }
        return result;
    }
    friend mpz_class fibonacci(const mpz_class &op);

    // void mpz_class::swap (mpz_class& op)
    // void swap (mpz_class& op1, mpz_class& op2)
    void swap(mpz_class &op) { mpz_swap(this->value, op.value); }
#if !defined ___GMPXX_STRICT_COMPATIBILITY___
    friend void swap(mpz_class &op1, mpz_class &op2) { mpz_swap(op1.value, op2.value); }
#endif
    friend std::ostream &operator<<(std::ostream &os, const mpz_class &op);
    friend std::ostream &operator<<(std::ostream &os, const mpz_t op);

    friend std::istream &operator>>(std::istream &stream, mpz_class &op);
    friend std::istream &operator>>(std::istream &stream, mpz_t op);

    // casts
    operator mpf_class() const;
    operator mpq_class() const;
    operator unsigned long int() const { return mpz_get_ui(this->value); }
    operator signed long int() const { return mpz_get_si(this->value); }
    operator unsigned int() const { return static_cast<unsigned int>(mpz_get_ui(this->value)); }
    operator signed int() const { return static_cast<signed int>(mpz_get_si(this->value)); }

    mpz_srcptr get_mpz_t() const { return value; }
    mpz_ptr get_mpz_t() { return value; }

  private:
    mpz_t value;
};
inline mpz_class &operator+=(mpz_class &op1, const mpz_class &op2) {
    mpz_add(op1.value, op1.value, op2.value);
    return op1;
}
inline mpz_class &operator-=(mpz_class &op1, const mpz_class &op2) {
    mpz_sub(op1.value, op1.value, op2.value);
    return op1;
}
inline mpz_class &operator*=(mpz_class &op1, const mpz_class &op2) {
    mpz_mul(op1.value, op1.value, op2.value);
    return op1;
}
inline mpz_class &operator/=(mpz_class &op1, const mpz_class &op2) {
    mpz_tdiv_q(op1.value, op1.value, op2.value);
    return op1;
}
inline mpz_class &operator&=(mpz_class &lhs, const mpz_class &rhs) {
    mpz_and(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpz_class &operator|=(mpz_class &lhs, const mpz_class &rhs) {
    mpz_ior(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpz_class &operator^=(mpz_class &lhs, const mpz_class &rhs) {
    mpz_xor(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpz_class &operator%=(mpz_class &op1, const mpz_class &op2) {
    mpz_tdiv_r(op1.value, op1.value, op2.value);
    return op1;
}
inline mpz_class operator+(const mpz_class &op) { return op; }
inline mpz_class operator-(const mpz_class &op) {
    mpz_class result;
    mpz_neg(result.value, op.value);
    return result;
}
inline mpz_class operator+(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_add(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator-(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_sub(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator*(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_mul(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator/(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_tdiv_q(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator%(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_tdiv_r(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator&(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_and(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator|(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_ior(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class operator^(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_xor(result.value, op1.value, op2.value);
    return result;
}

// +=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs) {
    mpz_add_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs) {
    if (rhs >= 0)
        mpz_add_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    else {
        mpz_sub_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
    }
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator+=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_add(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// +
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    mpz_add_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2) {
    mpz_class result(op2);
    mpz_add_ui(result.value, result.value, static_cast<unsigned long int>(op1));
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    if (op2 >= 0)
        mpz_add_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    else {
        mpz_sub_ui(result.value, result.value, static_cast<unsigned long int>(-op2));
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2) {
    mpz_class result(op2);
    if (op1 >= 0)
        mpz_add_ui(result.value, result.value, static_cast<unsigned long int>(op1));
    else
        mpz_sub_ui(result.value, result.value, static_cast<unsigned long int>(-op1));
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator+(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    result += op2;
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator+(const T op1, const mpz_class &op2) { return op2 + op1; }

// -=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs) {
    mpz_sub_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs) {
    if (rhs >= 0)
        mpz_sub_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    else {
        mpz_add_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
    }
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator-=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_sub(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// -
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2) {
    mpz_class result;
    mpz_sub_ui(result.value, op1.value, op2);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_ui_sub(result.value, op1, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    if (op2 >= 0)
        mpz_sub_ui(result.value, op1.value, static_cast<unsigned long int>(op2));
    else
        mpz_add_ui(result.value, op1.value, static_cast<unsigned long int>(-op2));

    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2) {
    mpz_class result;
    if (op1 >= 0) {
        mpz_ui_sub(result.value, static_cast<unsigned long int>(op1), op2.value);
    } else {
        mpz_add_ui(result.value, op2.value, static_cast<unsigned long int>(-op1));
        mpz_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator-(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_sub(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator-(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_sub(result.value, result.value, op2.value);
    return result;
}

// *=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs) {
    mpz_mul_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs) {
    mpz_mul_si(lhs.value, lhs.value, static_cast<signed long int>(rhs));
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator*=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_mul(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// *
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    mpz_mul_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2) {
    mpz_class result(op2);
    mpz_mul_ui(result.value, result.value, static_cast<unsigned long int>(op1));
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    mpz_mul_si(result.value, result.value, static_cast<signed long int>(op2));
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2) {
    mpz_class result(op2);
    mpz_mul_si(result.value, result.value, static_cast<signed long int>(op1));
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator*(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_mul(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator*(const T op1, const mpz_class &op2) { return op2 * op1; }

// /=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs) {
    mpz_tdiv_q_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs) {
    if (rhs >= 0)
        mpz_tdiv_q_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    else {
        mpz_tdiv_q_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
        mpz_neg(lhs.value, lhs.value);
    }
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator/=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_tdiv_q(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// /
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    mpz_tdiv_q_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_q(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    if (op2 >= 0)
        mpz_tdiv_q_ui(result.value, result.value, static_cast<signed long int>(op2));
    else {
        mpz_tdiv_q_ui(result.value, result.value, static_cast<unsigned long int>(-op2));
        mpz_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_q(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator/(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_tdiv_q(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator/(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_q(result.value, result.value, op2.value);
    return result;
}
// %=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs) {
    mpz_tdiv_r_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs) {
    if (rhs >= 0)
        mpz_tdiv_r_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    else {
        mpz_tdiv_r_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
    }
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator%=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_tdiv_r(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// %
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    mpz_tdiv_r_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_r(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2) {
    mpz_class result(op1);
    if (op2 >= 0)
        mpz_tdiv_r_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    else {
        mpz_tdiv_r_ui(result.value, result.value, static_cast<signed long int>(-op2));
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_r(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator%(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_tdiv_r(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator%(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_tdiv_r(result.value, result.value, op2.value);
    return result;
}
// &=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_and(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_and(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator&=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_and(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// &
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_and(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_and(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_and(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_and(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator&(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_and(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator&(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_and(result.value, result.value, op2.value);
    return result;
}
// |=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_ior(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_ior(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator|=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_ior(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// |
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_ior(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_ior(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_ior(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_ior(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator|(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_ior(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator|(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_ior(result.value, result.value, op2.value);
    return result;
}
// ^=
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_xor(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_xor(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpz_class &) operator^=(mpz_class &lhs, const T rhs) {
    mpz_class _rhs(rhs);
    mpz_xor(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
// ^
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_xor(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_xor(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_xor(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_xor(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator^(const mpz_class &op1, const T op2) {
    mpz_class result(op2);
    mpz_xor(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpz_class) operator^(const T op1, const mpz_class &op2) {
    mpz_class result(op1);
    mpz_xor(result.value, result.value, op2.value);
    return result;
}

/////
inline mpz_class &mpz_class::operator=(const signed long int op) {
    mpz_set_si(this->value, op);
    return *this;
}
inline mpz_class &mpz_class::operator=(const unsigned long int op) {
    mpz_set_ui(this->value, op);
    return *this;
}
inline mpz_class &mpz_class::operator=(const signed int op) {
    mpz_set_si(this->value, static_cast<signed long int>(op));
    return *this;
}
inline mpz_class &mpz_class::operator=(const unsigned int op) {
    mpz_set_ui(this->value, static_cast<unsigned long int>(op));
    return *this;
}
inline mpz_class &mpz_class::operator=(const signed char op) {
    mpz_set_si(this->value, static_cast<signed long int>(op));
    return *this;
}
inline mpz_class &mpz_class::operator=(const unsigned char op) {
    mpz_set_ui(this->value, static_cast<unsigned long int>(op));
    return *this;
}
inline mpz_class &mpz_class::operator=(const char op) {
    if (std::is_signed<char>::value)
        return *this = static_cast<signed char>(op);
    else
        return *this = static_cast<unsigned char>(op);
}
inline mpz_class abs(const mpz_class &op) {
    mpz_class result;
    mpz_abs(result.value, op.value);
    return result;
}
inline int sgn(const mpz_class &op) { return mpz_sgn(op.value); }
inline mpz_class sqrt(const mpz_class &op) {
    mpz_class result;
    mpz_sqrt(result.value, op.value);
    return result;
}
inline mpz_class gcd(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_gcd(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class lcm(const mpz_class &op1, const mpz_class &op2) {
    mpz_class result;
    mpz_lcm(result.value, op1.value, op2.value);
    return result;
}
inline mpz_class factorial(const mpz_class &n) {
    if (n < 0) {
        throw std::domain_error("factorial(negative)");
    }
    double log2_n = mpz_sizeinbase(n.get_mpz_t(), 2);
    if (log2_n > 300) {
        throw std::bad_alloc();
    }
    mpz_class result;
    try {
        mpz_fac_ui(result.get_mpz_t(), n.get_ui());
    } catch (const std::bad_alloc &) {
        throw;
    }
    return result;
}
inline mpz_class primorial(const mpz_class &op) {
    if (op < 0) {
        throw std::domain_error("primorial(negative)");
    }
    double log2_n = mpz_sizeinbase(op.get_mpz_t(), 2);
    if (log2_n > 300) {
        throw std::bad_alloc();
    }
    mpz_class result;
    try {
        mpz_primorial_ui(result.get_mpz_t(), op.get_ui());
    } catch (const std::bad_alloc &) {
        throw;
    }
    return result;
}
inline mpz_class fibonacci(const mpz_class &op) {
    double log2_op = mpz_sizeinbase(op.get_mpz_t(), 2);
    if (log2_op > 300) {
        throw std::bad_alloc();
    }
    mpz_class adjusted_op = op;
    bool isNegative = op < 0;
    if (isNegative) {
        adjusted_op = -op;
    }
    unsigned long int n = adjusted_op.get_ui();
    mpz_class result;
    mpz_fib_ui(result.get_mpz_t(), n);
    if (isNegative) {
        if ((op + 1) % 2 != 0) {
            result = -result;
        }
    }
    return result;
}
inline void print_mpz(std::ostream &os, const mpz_t op) {
    std::ios_base::fmtflags flags = os.flags();
    std::streamsize width = os.width();
    char fill = os.fill();
    char *str = nullptr;
    bool is_hex = flags & std::ios::hex;
    bool is_oct = flags & std::ios::oct;
    bool show_base = flags & std::ios::showbase;
    bool uppercase = flags & std::ios::uppercase;

    if (mpz_sgn(op) == 0) {
        if (is_hex && show_base) {
            str = strdup(uppercase ? "0X0" : "0x0");
        } else if (is_oct) {
            str = strdup("0");
        } else {
            str = strdup("0");
        }
    } else {
        if (is_oct) {
            gmp_asprintf(&str, (show_base) ? "%#Zo" : "%Zo", op);
        } else if (is_hex) {
            if (show_base && uppercase) {
                gmp_asprintf(&str, "%#ZX", op);
            } else if (show_base && !uppercase) {
                gmp_asprintf(&str, "%#Zx", op);
            } else if (!show_base && uppercase) {
                gmp_asprintf(&str, "%ZX", op);
            } else {
                gmp_asprintf(&str, "%Zx", op);
            }
        } else {
            gmp_asprintf(&str, "%Zd", op);
        }
    }
    std::string s(str);
    free(str);

    if (flags & std::ios::showpos && mpz_sgn(op) > 0) {
        s.insert(0, "+");
    }

    std::streamsize len = s.length();
    if (len < width) {
        std::streamsize padding_length = width - len;
        if (flags & std::ios::left) {
            s.append(padding_length, fill);
        } else if (flags & std::ios::internal) {
            size_t pos = 0;
            if (s[0] == '-' || s[0] == '+') {
                pos = 1;
            }
            if (s.length() > pos + 1 && (s[pos] == '0' && (s[pos + 1] == 'x' || s[pos + 1] == 'X'))) {
                pos += 2;
            }
            s.insert(pos, padding_length, fill);
        } else {
            s.insert(0, padding_length, fill);
        }
    }
    os << s;
    os.width(0);
}
inline std::ostream &operator<<(std::ostream &os, const mpz_class &op) {
    print_mpz(os, op.get_mpz_t());
    return os;
}
inline std::ostream &operator<<(std::ostream &os, const mpz_t &op) {
    print_mpz(os, op);
    return os;
}
inline std::istream &read_nofmtflags_mpz_from_stream(std::istream &stream, mpz_t op) {
    char ch;
    std::string number;
    bool negative = false;
    bool is_space = false;
    int base = 10;

    std::ios_base::fmtflags current_flags = stream.flags();
    int counter = 0;
    while (stream >> ch && isspace(ch)) {
        is_space = true;
        counter++;
    }
    if (!(current_flags & std::ios::skipws) && is_space == true) {
        for (int i = 0; i <= counter; i++)
            stream.unget();
        stream.setstate(std::ios::failbit);
        return stream;
    }
    if (ch == '+' || ch == '-') {
        negative = (ch == '-');
        stream.get(ch);
    }
    if (ch == '0') {
        stream.get(ch);
        if (ch == 'x' || ch == 'X') {
            base = 16;
            stream.get(ch);
        } else if (isdigit(ch) && ch != '8' && ch != '9') {
            base = 8;
        } else {
            // retake zero again
            stream.unget();
            stream.unget();
            stream.get(ch);
        }
    }
    while (isdigit(ch) || (base == 16 && isxdigit(ch))) {
        number += ch;
        if (!stream.get(ch))
            break;
    }
    if (!stream.eof())
        stream.unget();
    int ret = mpz_set_str(op, number.c_str(), base);
    if (ret != 0) {
        stream.setstate(stream.rdstate() & ~std::ios::goodbit);
        stream.setstate(std::ios::failbit);
    } else {
        stream.clear(stream.rdstate() & ~std::ios::failbit);
        stream.setstate(std::ios::goodbit);
    }
    if (negative) {
        mpz_neg(op, op);
    }
    return stream;
}
inline std::istream &read_base_mpz_from_stream(std::istream &stream, mpz_t op, int base = 10) {
    char ch;
    std::string number;
    bool negative = false;
    bool number_started = false;
    while (stream >> ch && isspace(ch))
        ;
    if (!stream) {
        stream.setstate(std::ios::failbit);
        return stream;
    }
    if (ch == '+' || ch == '-') {
        negative = (ch == '-');
        stream.get(ch);
    }
    while ((base <= 10 && isdigit(ch) && ch < '0' + base) || (base == 16 && isxdigit(ch))) {
        if (base == 16 && isalpha(ch)) {
            ch = tolower(ch);
        }
        number += ch;
        number_started = true;
        if (!stream.get(ch))
            break;
    }
    if (number_started && (ch < '0' || (ch > '9' && base <= 10) || (ch > 'f' && base == 16))) {
        stream.unget();
    }
    if (negative && !number.empty())
        number = "-" + number;
    if (number.empty()) {
        stream.setstate(std::ios::failbit);
    }
    int ret = mpz_set_str(op, number.c_str(), base);
    if (ret != 0) {
        stream.setstate(std::ios::failbit);
    } else {
        stream.clear(stream.rdstate() & ~std::ios::failbit);
        stream.setstate(std::ios::goodbit);
    }
    return stream;
}
inline std::istream &read_mpz_from_stream(std::istream &stream, mpz_t op) {
    std::ios_base::fmtflags current_flags = stream.flags();
    if (current_flags == std::ios_base::fmtflags(0)) {
        return read_nofmtflags_mpz_from_stream(stream, op);
    }
    if (current_flags & std::ios::oct) {
        return read_base_mpz_from_stream(stream, op, 8);
    } else if (current_flags & std::ios::hex) {
        return read_base_mpz_from_stream(stream, op, 16);
    }
    return read_base_mpz_from_stream(stream, op);
}
inline std::istream &operator>>(std::istream &stream, mpz_t op) { return read_mpz_from_stream(stream, op); }
inline std::istream &operator>>(std::istream &stream, mpz_class &op) { return read_mpz_from_stream(stream, op.get_mpz_t()); }
class mpq_class {
  public:
    ////////////////////////////////////////////////////////////////////////////////////////
    // 12.3 C++ Interface Rationals
    // cf. https://gmplib.org/manual/C_002b_002b-Interface-Rationals
    ////////////////////////////////////////////////////////////////////////////////////////
    // constructors and destructors
    mpq_class() { mpq_init(value); }
    // The rule of 0/3/5
    // The rule 1 of 5 copy constructor
    mpq_class(const mpq_class &op) {
        mpq_init(value);
        mpq_set(value, op.value);
    }
    // The rule 2 of 5 copy assignment operator
    mpq_class &operator=(const mpq_class &op) noexcept {
        if (this != &op) {
            mpq_set(value, op.value);
        }
        return *this;
    }
    // The rule 3 of 5 default deconstructor
    ~mpq_class() { mpq_clear(value); }
    // The rule 4 of 5 move constructor
    mpq_class(mpq_class &&op) noexcept {
        mpq_init(value);
        mpq_swap(value, op.value);
    }
    // The rule 5 of 5 move assignment operator
    mpq_class &operator=(mpq_class &&op) noexcept {
        if (this != &op) {
            mpq_swap(value, op.value);
        }
        return *this;
    }
    // constructors
    explicit mpq_class(const mpq_t q) {
        mpq_init(value);
        mpq_set(value, q);
    }
    mpq_class(const mpz_t op) {
        mpq_init(value);
        mpq_set_z(value, op);
    }
    mpq_class(const mpf_t op) {
        mpq_init(value);
        mpq_set_f(value, op);
    }
    mpq_class(const mpz_class &op1, const mpz_class &op2) {
        mpq_init(value);
        mpq_set_num(value, op1.get_mpz_t());
        mpq_set_den(value, op2.get_mpz_t());
        if (op2 == 0) {
            mpq_clear(value);
            throw std::invalid_argument("Denominator cannot be zero in a rational number.");
        }
        mpq_canonicalize(value);
    }
    mpq_class(const mpz_class &op) {
        mpq_init(value);
        mpq_set_z(value, op.get_mpz_t());
    }
    mpq_class(const char *str, int base = 0) {
        mpq_init(value);
        if (mpq_set_str(value, str, base) != 0) {
            throw std::invalid_argument("");
        }
    }
    mpq_class(const std::string &str, int base = 0) {
        mpq_init(value);
        if (mpq_set_str(value, str.c_str(), base) != 0) {
            throw std::invalid_argument("");
        }
    }
    mpq_class(unsigned long int op1, unsigned long int op2) {
        mpq_init(value);
        mpq_set_ui(value, op1, op2);
    }
    mpq_class(signed long int op1, signed long int op2) {
        mpq_init(value);
        mpq_set_si(value, op1, op2);
    }
    mpq_class(unsigned int op1, unsigned int op2) {
        mpq_init(value);
        mpq_set_ui(value, static_cast<unsigned long int>(op1), static_cast<unsigned long int>(op2));
    }
    mpq_class(int op1, int op2) {
        mpq_init(value);
        mpq_set_si(value, static_cast<signed long int>(op1), static_cast<signed long int>(op2));
    }
    mpq_class(unsigned long int op) {
        mpq_init(value);
        mpq_set_ui(value, op, static_cast<unsigned long int>(1));
    }
    mpq_class(signed long int op) {
        mpq_init(value);
        mpq_set_si(value, op, static_cast<unsigned long int>(1));
    }
    mpq_class(unsigned int op) {
        mpq_init(value);
        mpq_set_ui(value, static_cast<unsigned long int>(op), static_cast<unsigned long int>(1));
    }
    mpq_class(int op) {
        mpq_init(value);
        mpq_set_si(value, static_cast<signed long int>(op), static_cast<signed long int>(1));
    }
    mpq_class(double op) {
        mpq_init(value);
        mpq_set_d(value, op);
    }

    // assignments from other objects
    inline mpq_class &operator=(const mpz_class &op);
    inline mpq_class &operator=(const signed long int op);
    inline mpq_class &operator=(const unsigned long int op);
    inline mpq_class &operator=(const signed int op);
    inline mpq_class &operator=(const unsigned int op);
    inline mpq_class &operator=(const signed char op);
    inline mpq_class &operator=(const unsigned char op);
    inline mpq_class &operator=(const char op);
    inline mpq_class &operator=(const float op);
    inline mpq_class &operator=(const double op);
    inline mpq_class &operator=(const char *op);
    inline mpq_class &operator=(const std::string &op);

    // operators
    mpq_class &operator++() {
        mpq_class one = 1;
        mpq_add(value, value, one.value);
        return *this;
    }

    mpq_class operator++(int) {
        mpq_class original = *this;
        ++(*this);
        return original;
    }
    mpq_class &operator--() {
        mpq_class one = 1;
        mpq_sub(value, value, one.value);
        return *this;
    }
    mpq_class operator--(int) {
        mpq_class original = *this;
        --(*this);
        return original;
    }
    template <typename T> INT_COND(T, mpq_class &) operator<<=(T n) {
        mpq_mul_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> INT_COND(T, mpq_class &) operator>>=(T n) {
        mpq_div_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> friend INT_COND(T, mpq_class) operator<<(const mpq_class &op1, T op2) {
        mpq_class result(op1);
        mpq_mul_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }
    template <typename T> friend INT_COND(T, mpq_class) operator>>(const mpq_class &op1, T op2) {
        mpq_class result(op1);
        mpq_div_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }

    // mpq_class arithmetic operators
    inline friend mpq_class &operator+=(mpq_class &lhs, const mpq_class &rhs);
    inline friend mpq_class &operator-=(mpq_class &lhs, const mpq_class &rhs);
    inline friend mpq_class &operator*=(mpq_class &lhs, const mpq_class &rhs);
    inline friend mpq_class &operator/=(mpq_class &lhs, const mpq_class &rhs);
    inline friend mpq_class operator+(const mpq_class &op);
    inline friend mpq_class operator-(const mpq_class &op);
    inline friend mpq_class operator+(const mpq_class &op1, const mpq_class &op2);
    inline friend mpq_class operator-(const mpq_class &op1, const mpq_class &op2);
    inline friend mpq_class operator*(const mpq_class &op1, const mpq_class &op2);
    inline friend mpq_class operator/(const mpq_class &op1, const mpq_class &op2);

    // mpq_class comparison operators
    inline friend bool operator==(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) == 0; }
    inline friend bool operator!=(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) != 0; }
    inline friend bool operator<(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) < 0; }
    inline friend bool operator>(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) > 0; }
    inline friend bool operator<=(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) <= 0; }
    inline friend bool operator>=(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.value, op2.value) >= 0; }

    inline friend bool operator==(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) == 0; }
    inline friend bool operator!=(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) != 0; }
    inline friend bool operator<(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) < 0; }
    inline friend bool operator>(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) > 0; }
    inline friend bool operator<=(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) <= 0; }
    inline friend bool operator>=(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.value, op2.get_mpz_t()) >= 0; }
    inline friend bool operator==(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) == 0; }
    inline friend bool operator!=(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) != 0; }
    inline friend bool operator<(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) > 0; }
    inline friend bool operator>(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) < 0; }
    inline friend bool operator<=(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) >= 0; }
    inline friend bool operator>=(const mpz_class &op1, const mpq_class &op2) { return mpq_cmp_z(op2.value, op1.get_mpz_t()) <= 0; }

    // mpq_class comparison operators (template version)
    template <typename T> inline friend NON_GMP_COND(T, bool) operator==(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) == 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator!=(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) != 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator<(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) < 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator>(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) > 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator<=(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) <= 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator>=(const mpq_class &op1, T op2) { return mpq_cmp(op1.value, mpq_class(op2).get_mpq_t()) >= 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator==(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) == 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator!=(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) != 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator<(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) > 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator>(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) < 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator<=(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) >= 0; }
    template <typename T> inline friend NON_GMP_COND(T, bool) operator>=(T op1, const mpq_class &op2) { return mpq_cmp(op2.value, mpq_class(op1).get_mpq_t()) <= 0; }

    // mpq_class arithmetic operators (template version)
    template <typename T> inline friend mpq_class &operator+=(mpq_class &lhs, const T rhs) {
        mpq_class _rhs(rhs);
        mpq_add(lhs.value, lhs.value, _rhs.value);
        return lhs;
    }
    template <typename T> inline friend mpq_class operator+(const mpq_class &op1, const T op2) {
        mpq_class result(op2);
        mpq_add(result.value, op1.value, result.value);
        return result;
    }
    template <typename T> inline friend mpq_class operator+(const T op1, const mpq_class &op2) { return op2 + op1; }
    template <typename T> inline friend mpq_class &operator-=(mpq_class &lhs, const T rhs) {
        mpq_class _rhs(rhs);
        mpq_sub(lhs.value, lhs.value, _rhs.value);
        return lhs;
    }
    template <typename T> inline friend mpq_class operator-(const mpq_class &op1, const T op2) {
        mpq_class result(op2);
        mpq_sub(result.value, op1.value, result.value);
        return result;
    }
    template <typename T> inline friend mpq_class operator-(const T op1, const mpq_class &op2) {
        mpq_class result(op1);
        mpq_sub(result.value, result.value, op2.value);
        return result;
    }
    template <typename T> inline friend mpq_class &operator*=(mpq_class &lhs, const T rhs) {
        mpq_class _rhs(rhs);
        mpq_mul(lhs.value, lhs.value, _rhs.value);
        return lhs;
    }
    template <typename T> inline friend mpq_class operator*(const mpq_class &op1, const T op2) {
        mpq_class result(op2);
        mpq_mul(result.value, op1.value, result.value);
        return result;
    }
    template <typename T> inline friend mpq_class operator*(const T op1, const mpq_class &op2) { return op2 * op1; }
    template <typename T> inline friend mpq_class &operator/=(mpq_class &lhs, const T rhs) {
        mpq_class _rhs(rhs);
        mpq_div(lhs.value, lhs.value, _rhs.value);
        return lhs;
    }
    template <typename T> inline friend mpq_class operator/(const mpq_class &op1, const T op2) {
        mpq_class result(op2);
        mpq_div(result.value, op1.value, result.value);
        return result;
    }
    template <typename T> inline friend mpq_class operator/(const T op1, const mpq_class &op2) {
        mpq_class result(op1);
        mpq_div(result.value, result.value, op2.value);
        return result;
    }
    // void mpq_class::canonicalize ()
    // mpq_class abs (mpq_class op)
    // double mpq_class::get_d (void)
    void canonicalize() { mpq_canonicalize(value); }
    friend mpq_class abs(const mpq_class &op);
    double get_d() const { return mpq_get_d(value); }

    // string mpq_class::get_str (int base = 10)
    // int mpq_class::set_str (const char *str, int base)
    // int mpq_class::set_str (const string& str, int base)
    std::string get_str(int base = 10) const {
        char *str = mpq_get_str(NULL, base, value);
        std::string result(str);
        void (*freefunc)(void *, size_t);
        mp_get_memory_functions(NULL, NULL, &freefunc);
        freefunc(str, std::strlen(str) + 1);
        return result;
    }
    int set_str(const char *str, int base = 10) {
        int ret = mpq_set_str(value, str, base);
        if (ret == 0) {
            mpq_canonicalize(value);
        }
        return ret;
    }
    int set_str(const std::string &str, int base = 10) { return set_str(str.c_str(), base); }
    // int sgn (mpq_class op)
    // void mpq_class::swap (mpq_class& op)
    // void swap (mpq_class& op1, mpq_class& op2)
    void swap(mpq_class &op) { mpq_swap(this->value, op.value); }
    friend int sgn(const mpq_class &op) { return mpq_sgn(op.value); }
#if !defined ___GMPXX_STRICT_COMPATIBILITY___
    friend void swap(mpq_class &op1, mpq_class &op2) { mpq_swap(op1.value, op2.value); }
#endif

    // mpz_class& mpq_class::get_num ()
    // mpz_class& mpq_class::get_den ()
    // mpz_t mpq_class::get_num_mpz_t ()
    // mpz_t mpq_class::get_den_mpz_t ()
    mpz_class get_num() const {
        mpz_class num(mpq_numref(value));
        return num;
    }
    mpz_class get_den() const {
        mpz_class den(mpq_denref(value));
        return den;
    }
    mpz_srcptr get_num_mpz_t() const { return mpq_numref(value); }
    mpz_srcptr get_den_mpz_t() const { return mpq_denref(value); }
    mpz_ptr get_num_mpz_t() { return mpq_numref(value); }
    mpz_ptr get_den_mpz_t() { return mpq_denref(value); }

    // istream& operator>> (istream& stream, mpq_class& rop)
    friend std::ostream &operator<<(std::ostream &os, const mpq_class &op);
    friend std::ostream &operator<<(std::ostream &os, const mpq_t op);

    friend std::istream &operator>>(std::istream &stream, mpq_class &op);
    friend std::istream &operator>>(std::istream &stream, mpq_t op);

    operator mpf_class() const;
    operator mpz_class() const;
    mpq_srcptr get_mpq_t() const { return value; }
    mpq_ptr get_mpq_t() { return value; }

  private:
    mpq_t value;
};

inline mpq_class &mpq_class::operator=(const mpz_class &op) {
    mpq_set_z(this->value, op.get_mpz_t());
    return *this;
}
inline mpq_class &mpq_class::operator=(signed long int op) {
    mpq_set_si(this->value, op, (signed long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(unsigned long int op) {
    mpq_set_ui(this->value, op, (unsigned long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(signed int op) {
    mpq_set_si(this->value, (signed long int)op, (signed long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(unsigned int op) {
    mpq_set_ui(this->value, (unsigned long int)op, (unsigned long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(signed char op) {
    mpq_set_si(this->value, (signed long int)op, (signed long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(unsigned char op) {
    mpq_set_ui(this->value, (unsigned long int)op, (unsigned long int)1);
    return *this;
}
inline mpq_class &mpq_class::operator=(char op) {
    if (std::is_signed<char>::value)
        return *this = static_cast<signed char>(op);
    else
        return *this = static_cast<unsigned char>(op);
}
inline mpq_class &mpq_class::operator=(float op) {
    mpq_set_d(this->value, (double)op);
    return *this;
}
inline mpq_class &mpq_class::operator=(double op) {
    mpq_set_d(this->value, op);
    return *this;
}
inline mpq_class &mpq_class::operator=(const char *op) {
    if (mpq_set_str(value, op, 10) != 0) {
        throw std::invalid_argument("Invalid string format for mpq_class");
    }
    return *this;
}
inline mpq_class &mpq_class::operator=(const std::string &op) {
    if (mpq_set_str(value, op.c_str(), 10) != 0) {
        throw std::invalid_argument("Invalid string format for mpq_class");
    }
    return *this;
}

inline mpq_class &operator+=(mpq_class &op1, const mpq_class &op2) {
    mpq_add(op1.value, op1.value, op2.value);
    return op1;
}
inline mpq_class &operator-=(mpq_class &op1, const mpq_class &op2) {
    mpq_sub(op1.value, op1.value, op2.value);
    return op1;
}
inline mpq_class &operator/=(mpq_class &op1, const mpq_class &op2) {
    mpq_div(op1.value, op1.value, op2.value);
    return op1;
}
inline mpq_class &operator*=(mpq_class &op1, const mpq_class &op2) {
    mpq_mul(op1.value, op1.value, op2.value);
    return op1;
}
inline mpq_class operator+(const mpq_class &op) { return op; }
inline mpq_class operator-(const mpq_class &op) {
    mpq_class result;
    mpq_neg(result.value, op.value);
    return result;
}
inline mpq_class operator+(const mpq_class &op1, const mpq_class &op2) {
    mpq_class result;
    mpq_add(result.value, op1.value, op2.value);
    return result;
}
inline mpq_class operator-(const mpq_class &op1, const mpq_class &op2) {
    mpq_class result;
    mpq_sub(result.value, op1.value, op2.value);
    return result;
}
inline mpq_class operator*(const mpq_class &op1, const mpq_class &op2) {
    mpq_class result;
    mpq_mul(result.value, op1.value, op2.value);
    return result;
}
inline mpq_class operator/(const mpq_class &op1, const mpq_class &op2) {
    mpq_class result;
    mpq_div(result.value, op1.value, op2.value);
    return result;
}
inline mpq_class abs(const mpq_class &op) {
    mpq_class rop(op);
    mpq_abs(rop.value, op.get_mpq_t());
    return rop;
}
inline void print_mpq(std::ostream &os, const mpq_t op) {
    std::ios_base::fmtflags flags = os.flags();
    std::streamsize width = os.width();
    char fill = os.fill();
    char *str = nullptr;
    bool is_hex = flags & std::ios::hex;
    bool is_oct = flags & std::ios::oct;
    bool show_base = flags & std::ios::showbase;
    bool uppercase = flags & std::ios::uppercase;
    std::string format;

    mpz_class den(mpq_denref(op));
    mpz_class num(mpq_numref(op));

    if ((num == 0 && den == 1) || (num == 0 && den == 0)) {
        if (is_oct) {         // is_oct, (show_base can be ignored since octal 0 = 0 or 00, and we use 0).
            if (width == 0) { // is_oct, width==0
                str = strdup("0");
            } else { // is_oct, width!=0
                str = strdup("0/0");
            }
        } else if (is_hex) {
            if (show_base) {      // is_hex, show_base
                if (width == 0) { // is_hex, show_base, width==0
                    str = strdup(uppercase ? "0X0" : "0x0");
                } else {
                    if (uppercase) {
                        str = strdup("0X0/0X0");
                    } else {
                        str = strdup("0x0/0x0");
                    }
                }
            } else {              // is_hex
                if (width == 0) { // is_hex, width==0
                    str = strdup("0");
                } else { // is_hex, width!=0
                    str = strdup("0/0");
                }
            }
        } else {              // is_dec
            if (width == 0) { // is_dec, width==0
                str = strdup("0");
            } else { // is_dec, width!=0
                str = strdup("0/0");
            }
        }
    } else if (den == 0) {
        if (is_oct) {        // is_oct, (show_base can be ignored since octal 0 = 0 or 00, and we use 0).
            if (show_base) { // is_oct, show_base
                gmp_asprintf(&str, "%#Qo", op);
            } else { // is_oct
                gmp_asprintf(&str, "%Qo", op);
            }
        } else if (is_hex) {
            if (show_base) { // is_hex, show_base
                if (uppercase) {
                    gmp_asprintf(&str, "%#QX", op);
                } else {
                    gmp_asprintf(&str, "%#Qx", op);
                }
            } else { // is_hex
                if (uppercase) {
                    gmp_asprintf(&str, "%QX", op);
                } else {
                    gmp_asprintf(&str, "%Qx", op);
                }
            }
            // Add 'x0' to "/0" to make it "/0x0"
            char *slashZero = strstr(str, "/0");
            size_t newLen = strlen(str) + 2;
            char *newStr = (char *)malloc(newLen + 1);
            if (!newStr) {
                free(str);
                throw std::bad_alloc();
            }
            size_t offset = slashZero - str;
            strncpy(newStr, str, offset);
            newStr[offset] = '\0';

            strcat(newStr, uppercase ? "/0X0" : "/0x0");
            strcat(newStr, slashZero + 2);
            free(str);
            str = newStr;
        } else { // is_dec
            gmp_asprintf(&str, "%Qd", op);
        }
    } else if (num == 0 && den != 1 && den != 0) {
        if (is_oct) {        // is_oct, (show_base can be ignored since octal 0 = 0 or 00, and we use 0).
            if (show_base) { // is_oct, show_base
                gmp_asprintf(&str, "%#Qo", op);
            } else { // is_oct
                gmp_asprintf(&str, "%Qo", op);
            }
        } else if (is_hex) {
            if (show_base) { // is_hex, show_base
                if (uppercase) {
                    gmp_asprintf(&str, "%#QX", op);
                } else {
                    gmp_asprintf(&str, "%#Qx", op);
                }
                // Add 'x0' to "0/" to make it "0x0/"
                char *zeroSlash = strstr(str, "0/");
                if (zeroSlash) {
                    size_t newLen = strlen(str) + 2;
                    char *newStr = (char *)malloc(newLen + 1);
                    if (!newStr) {
                        free(str);
                        throw std::bad_alloc();
                    }
                    size_t offset = zeroSlash - str;
                    strncpy(newStr, str, offset + 1);
                    newStr[offset + 1] = '\0';
                    strcat(newStr, uppercase ? "X0/" : "x0/");
                    strcat(newStr, zeroSlash + 2);
                    free(str);
                    str = newStr;
                }
            } else { // is_hex
                if (uppercase) {
                    gmp_asprintf(&str, "%QX", op);
                } else {
                    gmp_asprintf(&str, "%Qx", op);
                }
            }
        } else { // is_dec
            gmp_asprintf(&str, "%Qd", op);
        }
    } else {
        if (is_oct) {        // is_oct, (show_base can be ignored since octal 0 = 0 or 00, and we use 0).
            if (show_base) { // is_oct, show_base
                gmp_asprintf(&str, "%#Qo", op);
            } else { // is_oct
                gmp_asprintf(&str, "%Qo", op);
            }
        } else if (is_hex) {
            if (show_base) { // is_hex, show_base
                if (uppercase) {
                    gmp_asprintf(&str, "%#QX", op);
                } else {
                    gmp_asprintf(&str, "%#Qx", op);
                }
            } else { // is_hex
                if (uppercase) {
                    gmp_asprintf(&str, "%QX", op);
                } else {
                    gmp_asprintf(&str, "%Qx", op);
                }
            }
        } else { // is_dec
            gmp_asprintf(&str, "%Qd", op);
        }
    }
    std::string s(str);
    free(str);

    if (flags & std::ios::showpos && mpq_sgn(op) > 0) {
        s.insert(0, "+");
    }
    std::streamsize len = s.length();
    if (len < width) {
        std::streamsize padding_length = width - len;
        if (flags & std::ios::left) {
            s.append(padding_length, fill);
        } else if (flags & std::ios::internal) {
            size_t pos = 0;
            if (s[0] == '-' || s[0] == '+') {
                pos = 1;
            }
            if (s.length() > pos + 1 && (s[pos] == '0' && (s[pos + 1] == 'x' || s[pos + 1] == 'X'))) {
                pos += 2;
            }
            s.insert(pos, padding_length, fill);
        } else {
            s.insert(0, padding_length, fill);
        }
    }
    os << s;
    os.width(0);
}
inline std::ostream &operator<<(std::ostream &os, const mpq_class &op) {
    print_mpq(os, op.get_mpq_t());
    return os;
}
inline std::ostream &operator<<(std::ostream &os, const mpq_t &op) {
    print_mpq(os, op);
    return os;
}
inline std::istream &read_base_mpq_from_stream(std::istream &stream, mpq_t op, int base = 10) {
    char ch;
    std::string number;
    mpz_ptr _numerator = mpq_numref(op);
    mpz_ptr _denominator = mpq_denref(op);
    mpq_class result;
    bool has_slash = false;

    read_base_mpz_from_stream(stream, _numerator, base);
    if (!stream.eof()) {
        while (true) {
            if (stream.get(ch)) {
                if (ch == '/') {
                    has_slash = true;
                    break;
                } else {
                    stream.unget();
                    break;
                }
            } else {
                break;
            }
        }
    }
    if (has_slash) {
        char _ch = stream.peek();
        if (isxdigit(_ch)) {
            read_base_mpz_from_stream(stream, _denominator, base);
        } else {
            stream.setstate(std::ios::failbit);
        }
    } else {
        mpz_set_ui(_denominator, 1UL);
    }
    return stream;
}
inline std::istream &read_nofmtflags_mpq_from_stream(std::istream &stream, mpq_t op) {
    char ch;
    std::string number;
    mpz_ptr _numerator = mpq_numref(op);
    mpz_ptr _denominator = mpq_denref(op);
    mpq_class result;
    bool has_slash = false;
    read_nofmtflags_mpz_from_stream(stream, _numerator);
    if (!stream.eof()) {
        stream.get(ch);
        if (ch == '/') {
            has_slash = true;
        } else {
            stream.unget();
        }
    }
    if (has_slash) {
        char _ch = stream.peek();
        if (isxdigit(_ch)) {
            read_nofmtflags_mpz_from_stream(stream, _denominator);
        } else {
            stream.setstate(std::ios::failbit);
        }
    } else
        mpz_set_ui(_denominator, 1UL);
    return stream;
}
inline std::istream &read_mpq_from_stream(std::istream &stream, mpq_t op) {
    std::ios_base::fmtflags current_flags = stream.flags();
    if (current_flags == std::ios_base::fmtflags(0)) {
        return read_nofmtflags_mpq_from_stream(stream, op);
    }
    if (current_flags & std::ios::oct) {
        return read_base_mpq_from_stream(stream, op, 8);
    } else if (current_flags & std::ios::hex) {
        return read_base_mpq_from_stream(stream, op, 16);
    }
    return read_base_mpq_from_stream(stream, op);
}
inline std::istream &operator>>(std::istream &stream, mpq_t op) { return read_mpq_from_stream(stream, op); }
inline std::istream &operator>>(std::istream &stream, mpq_class &op) { return read_mpq_from_stream(stream, op.get_mpq_t()); }
class mpf_class {
  public:
    ////////////////////////////////////////////////////////////////////////////////////////
    // 12.4 C++ Interface Floats
    // https://gmplib.org/manual/C_002b_002b-Interface-Floats
    ////////////////////////////////////////////////////////////////////////////////////////
    // constructors and destructors
    mpf_class() { mpf_init(value); } // a default precision has already been established during initialization by a call to mpf_set_default_prec.
    // The rule of 0/3/5
    // The rule 1 of 5 copy constructor
    mpf_class(const mpf_class &op) {
        mpf_init2(value, mpf_get_prec(op.value));
        mpf_set(value, op.value);
    }
    mpf_class(const mpf_class &op, mp_bitcnt_t prec) {
        mpf_init2(value, prec);
        mpf_set(value, op.value);
    }
    // The rule 2 of 5 copy assignment operator
    mpf_class &operator=(const mpf_class &op) noexcept {
        if (this != &op) {
            mpf_set(value, op.value);
        }
        return *this;
    }
    // The rule 3 of 5 default deconstructor
    ~mpf_class() { mpf_clear(value); }
    // The rule 4 of 5 move constructor
    mpf_class(mpf_class &&op) noexcept {
        mpf_init(value);
        mpf_swap(value, op.value);
    }
    // The rule 5 of 5 move assignment operator
    mpf_class &operator=(mpf_class &&op) noexcept {
        if (this != &op) {
#if !defined ___GMPXX_MKII_NOPRECCHANGE___
            if (mpf_get_prec(this->get_mpf_t()) == mpf_get_prec(op.value)) {
                mpf_swap(value, op.value);
            } else {
                mpf_set(value, op.value);
            }
#else
            mpf_swap(value, op.value);
#endif
        }
        return *this;
    }
    // constructors
    explicit mpf_class(const mpf_t op) {
        mp_bitcnt_t op_prec = mpf_get_prec(op);
        mpf_init2(value, op_prec);
        mpf_set(value, op);
    }
    mpf_class(const mpz_t op) noexcept {
        mpf_init(value);
        mpf_set_z(value, op);
    }
    mpf_class(const mpq_t op) noexcept {
        mpf_init(value);
        mpf_set_q(value, op);
    }
    mpf_class(const unsigned long int op) noexcept { mpf_init_set_ui(value, op); }
    mpf_class(const unsigned int op) noexcept : mpf_class(static_cast<unsigned long int>(op)) {}
    mpf_class(const signed long int op) noexcept { mpf_init_set_si(value, op); }
    mpf_class(const signed int op) noexcept { mpf_init_set_si(value, static_cast<signed long int>(op)); }
    mpf_class(const double op) noexcept { mpf_init_set_d(value, op); }
    mpf_class(const char *str) {
        if (mpf_init_set_str(value, str, gmpxx_defaults::base) != 0) {
            mpf_clear(value);
            throw std::invalid_argument("");
        }
    }
    mpf_class(const std::string &str) {
        if (mpf_init_set_str(value, str.c_str(), gmpxx_defaults::base) != 0) {
            mpf_clear(value);
            throw std::invalid_argument("");
        }
    }

    explicit mpf_class(const mpf_t op, mp_bitcnt_t prec) {
        mpf_init2(value, prec);
        mpf_set(value, op);
    }
    mpf_class(const mpz_t op, mp_bitcnt_t prec) {
        mpf_init2(value, prec);
        mpf_set_z(value, op);
    }
    mpf_class(const mpq_t op, mp_bitcnt_t prec) {
        mpf_init2(value, prec);
        mpf_set_q(value, op);
    }
    mpf_class(const unsigned long int op, mp_bitcnt_t prec) noexcept {
        mpf_init2(value, prec);
        mpf_set_ui(value, op);
    }
    mpf_class(const unsigned int op, mp_bitcnt_t prec) noexcept {
        mpf_init2(value, prec);
        mpf_set_ui(value, static_cast<unsigned long int>(op));
    }
    mpf_class(const signed long int op, mp_bitcnt_t prec) noexcept {
        mpf_init2(value, prec);
        mpf_set_si(value, op);
    }
    mpf_class(const signed int op, mp_bitcnt_t prec) noexcept {
        mpf_init2(value, prec);
        mpf_set_si(value, static_cast<signed long int>(op));
    }
    mpf_class(const double op, mp_bitcnt_t prec) noexcept {
        mpf_init2(value, prec);
        mpf_set_d(value, op);
    }
    mpf_class(const char *str, mp_bitcnt_t prec, int base = gmpxx_defaults::base) {
        mpf_init2(value, prec);
        if (mpf_set_str(value, str, base) != 0) {
            throw std::invalid_argument("");
        }
    }
    mpf_class(const std::string &str, mp_bitcnt_t prec, int base = gmpxx_defaults::base) {
        mpf_init2(value, prec);
        if (mpf_set_str(value, str.c_str(), base) != 0) {
            throw std::invalid_argument("");
        }
    }

    // assignments from other objects
    mpf_class &operator=(double d) noexcept {
        mpf_set_d(value, d);
        return *this;
    }
    mpf_class &operator=(unsigned long int d) noexcept {
        mpf_set_ui(value, d);
        return *this;
    }
    mpf_class &operator=(signed long int d) noexcept {
        mpf_set_si(value, d);
        return *this;
    }
    mpf_class &operator=(unsigned int d) noexcept {
        mpf_set_ui(value, static_cast<unsigned long int>(d));
        return *this;
    }
    mpf_class &operator=(signed int d) noexcept {
        mpf_set_si(value, static_cast<signed long int>(d));
        return *this;
    }
    mpf_class &operator=(const char *str) {
        if (mpf_set_str(value, str, gmpxx_defaults::base) != 0) {
            throw std::invalid_argument("");
        }
        return *this;
    }
    mpf_class &operator=(const std::string &str) {
        if (mpf_set_str(value, str.c_str(), gmpxx_defaults::base) != 0) {
            throw std::invalid_argument("");
        }
        return *this;
    }
    // operators
    inline mpf_class &operator++() {
        mpf_add_ui(value, value, 1);
        return *this;
    }
    inline mpf_class &operator--() {
        mpf_sub_ui(value, value, 1);
        return *this;
    }
    inline mpf_class operator++(int) {
        mpf_add_ui(value, value, 1);
        return *this;
    }
    inline mpf_class operator--(int) {
        mpf_sub_ui(value, value, 1);
        return *this;
    }
    template <typename T> INT_COND(T, mpf_class &) operator<<=(T n) {
        mpf_mul_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> INT_COND(T, mpf_class &) operator>>=(T n) {
        mpf_div_2exp(value, value, static_cast<mp_bitcnt_t>(n));
        return *this;
    }
    template <typename T> friend INT_COND(T, mpf_class) operator<<(const mpf_class &op1, T op2) {
        mpf_class result(op1);
        mpf_mul_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }
    template <typename T> friend INT_COND(T, mpf_class) operator>>(const mpf_class &op1, T op2) {
        mpf_class result(op1);
        mpf_div_2exp(result.value, result.value, static_cast<mp_bitcnt_t>(op2));
        return result;
    }

    // mpf_class arithmetic operators
    inline friend mpf_class &operator+=(mpf_class &lhs, const mpf_class &rhs);
    inline friend mpf_class &operator-=(mpf_class &lhs, const mpf_class &rhs);
    inline friend mpf_class &operator*=(mpf_class &lhs, const mpf_class &rhs);
    inline friend mpf_class &operator/=(mpf_class &lhs, const mpf_class &rhs);
    inline friend mpf_class operator+(const mpf_class &op);
    inline friend mpf_class operator-(const mpf_class &op);
    inline friend mpf_class operator+(const mpf_class &op1, const mpf_class &op2);
    inline friend mpf_class operator-(const mpf_class &op1, const mpf_class &op2);
    inline friend mpf_class operator*(const mpf_class &op1, const mpf_class &op2);
    inline friend mpf_class operator/(const mpf_class &op1, const mpf_class &op2);

    inline friend mpf_class &operator+=(mpf_class &lhs, const mpz_class &rhs);
    inline friend mpf_class &operator-=(mpf_class &lhs, const mpz_class &rhs);
    inline friend mpf_class &operator*=(mpf_class &lhs, const mpz_class &rhs);
    inline friend mpf_class &operator/=(mpf_class &lhs, const mpz_class &rhs);
    inline friend mpf_class operator+(const mpf_class &op1, const mpz_class &op2);
    inline friend mpf_class operator+(const mpz_class &op1, const mpf_class &op2);
    inline friend mpf_class operator-(const mpf_class &op1, const mpz_class &op2);
    inline friend mpf_class operator-(const mpz_class &op1, const mpf_class &op2);
    inline friend mpf_class operator*(const mpf_class &op1, const mpz_class &op2);
    inline friend mpf_class operator*(const mpz_class &op1, const mpf_class &op2);
    inline friend mpf_class operator/(const mpf_class &op1, const mpz_class &op2);
    inline friend mpf_class operator/(const mpz_class &op1, const mpf_class &op2);

    inline friend mpf_class &operator+=(mpf_class &lhs, const mpq_class &rhs);
    inline friend mpf_class &operator-=(mpf_class &lhs, const mpq_class &rhs);
    inline friend mpf_class &operator*=(mpf_class &lhs, const mpq_class &rhs);
    inline friend mpf_class &operator/=(mpf_class &lhs, const mpq_class &rhs);
    inline friend mpf_class operator+(const mpf_class &op1, const mpq_class &op2);
    inline friend mpf_class operator+(const mpq_class &op1, const mpf_class &op2);
    inline friend mpf_class operator-(const mpf_class &op1, const mpq_class &op2);
    inline friend mpf_class operator-(const mpq_class &op1, const mpf_class &op2);
    inline friend mpf_class operator*(const mpf_class &op1, const mpq_class &op2);
    inline friend mpf_class operator*(const mpq_class &op1, const mpf_class &op2);
    inline friend mpf_class operator/(const mpf_class &op1, const mpq_class &op2);
    inline friend mpf_class operator/(const mpq_class &op1, const mpf_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2);

    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend NON_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2);
    template <typename T> inline friend UNSIGNED_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2);
    template <typename T> inline friend SIGNED_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2);
    template <typename T> inline friend NON_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2);

    // mpf_class comparison operators
    inline friend bool operator==(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) == 0; }
    inline friend bool operator!=(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) != 0; }
    inline friend bool operator<(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) < 0; }
    inline friend bool operator>(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) > 0; }
    inline friend bool operator<=(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) <= 0; }
    inline friend bool operator>=(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.value, op2.value) >= 0; }

    inline friend bool operator==(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) == 0; }
    inline friend bool operator!=(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) != 0; }
    inline friend bool operator<(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) < 0; }
    inline friend bool operator>(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) > 0; }
    inline friend bool operator<=(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) <= 0; }
    inline friend bool operator>=(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.value, op2.get_mpz_t()) >= 0; }

    inline friend bool operator==(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) == 0; }
    inline friend bool operator!=(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) != 0; }
    inline friend bool operator<(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) > 0; }
    inline friend bool operator>(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) < 0; }
    inline friend bool operator<=(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) >= 0; }
    inline friend bool operator>=(const mpz_class &op1, const mpf_class &op2) { return mpf_cmp_z(op2.value, op1.get_mpz_t()) <= 0; }

    inline friend bool operator==(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) == 0; }
    inline friend bool operator!=(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) != 0; }
    inline friend bool operator<(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) < 0; }
    inline friend bool operator>(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) > 0; }
    inline friend bool operator<=(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) <= 0; }
    inline friend bool operator>=(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.value, mpf_class(op2).value) >= 0; }

    inline friend bool operator==(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) == 0; }
    inline friend bool operator!=(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) != 0; }
    inline friend bool operator<(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) > 0; }
    inline friend bool operator>(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) < 0; }
    inline friend bool operator<=(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) >= 0; }
    inline friend bool operator>=(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(op2.value, mpf_class(op1).value) <= 0; }

    inline friend bool operator==(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) == 0; }
    inline friend bool operator!=(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) != 0; }
    inline friend bool operator<(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) < 0; }
    inline friend bool operator>(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) > 0; }
    inline friend bool operator<=(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) <= 0; }
    inline friend bool operator>=(const mpf_class &op1, double op2) { return mpf_cmp_d(op1.value, op2) >= 0; }

    inline friend bool operator==(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) == 0; }
    inline friend bool operator!=(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) != 0; }
    inline friend bool operator<(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) > 0; }
    inline friend bool operator>(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) < 0; }
    inline friend bool operator<=(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) >= 0; }
    inline friend bool operator>=(double op1, const mpf_class &op2) { return mpf_cmp_d(op2.value, op1) <= 0; }

    // mpf_class comparison operators (template version)
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator==(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) == 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator!=(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) != 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) < 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) > 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<=(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) <= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>=(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.value, static_cast<unsigned long int>(op2)) >= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator==(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) == 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator!=(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) != 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) > 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) < 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator<=(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) >= 0; }
    template <typename T> inline friend UNSIGNED_INT_COND(T, bool) operator>=(T op1, const mpf_class &op2) { return mpf_cmp_ui(op2.value, static_cast<unsigned long int>(op1)) <= 0; }

    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator==(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) == 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator!=(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) != 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) < 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) > 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<=(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) <= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>=(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.value, static_cast<signed long int>(op2)) >= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator==(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) == 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator!=(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) != 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) > 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) < 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator<=(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) >= 0; }
    template <typename T> inline friend SIGNED_INT_COND(T, bool) operator>=(T op1, const mpf_class &op2) { return mpf_cmp_si(op2.value, static_cast<signed long int>(op1)) <= 0; }

    // mpf_class abs (mpf_class op)
    // mpf_class ceil (mpf_class op)
    friend mpf_class abs(const mpf_class &op);
    friend mpf_class ceil(const mpf_class &op);

    // bool mpf_class::fits_sint_p (void)
    // bool mpf_class::fits_slong_p (void)
    // bool mpf_class::fits_sshort_p (void)
    // bool mpf_class::fits_uint_p (void)
    // bool mpf_class::fits_ulong_p (void)
    // bool mpf_class::fits_ushort_p (void)
    bool fits_sint_p() const { return mpf_fits_sint_p(value); }
    bool fits_slong_p() const { return mpf_fits_slong_p(value); }
    bool fits_sshort_p() const { return mpf_fits_sshort_p(value); }
    bool fits_uint_p() const { return mpf_fits_uint_p(value); }
    bool fits_ulong_p() const { return mpf_fits_ulong_p(value); }
    bool fits_ushort_p() const { return mpf_fits_ushort_p(value); }

    // mpf_class floor (mpf_class op)
    // mpf_class hypot (mpf_class op1, mpf_class op2)
    friend mpf_class floor(const mpf_class &op);
    friend mpf_class hypot(const mpf_class &op1, const mpf_class &op2);

    // double mpf_class::get_d (void)
    // long mpf_class::get_si (void)
    // unsigned long mpf_class::get_ui (void)
    // string mpf_class::get_str (mp_exp_t& exp, int base = 10, size_t digits = 0)
    double get_d() const { return mpf_get_d(value); }
    unsigned long get_ui() const { return mpf_get_ui(value); }
    long get_si() const { return mpf_get_si(value); }
    std::string get_str(mp_exp_t &exp, int base = 10, size_t digits = 0) const {
        char *temp = mpf_get_str(nullptr, &exp, base, digits, value);
        std::string result = temp;
        void (*freefunc)(void *, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(temp, std::strlen(temp) + 1);
        return result;
    }
    void div_2exp(mp_bitcnt_t exp) {
        mpf_ptr non_const_ptr = const_cast<mpf_ptr>(this->get_mpf_t());
        mpf_div_2exp(non_const_ptr, this->get_mpf_t(), exp);
    }
    void mul_2exp(mp_bitcnt_t exp) {
        mpf_ptr non_const_ptr = const_cast<mpf_ptr>(this->get_mpf_t());
        mpf_mul_2exp(non_const_ptr, this->get_mpf_t(), exp);
    }
    void set_epsilon() {
        mp_bitcnt_t bits = get_prec();
        mpf_set_ui(this->get_mpf_t(), 1);                             // epsilon = 1
        mpf_div_2exp(this->get_mpf_t(), this->get_mpf_t(), bits - 1); // epsilon = 1 / 2^(bits-1)
    }
    // int mpf_class::set_str (const char *str, int base)
    // int mpf_class::set_str (const string& str, int base)
    int set_str(const char *str, int base) { return mpf_set_str(value, str, base); }
    int set_str(const std::string &str, int base) { return mpf_set_str(value, str.c_str(), base); }

    // int sgn (mpf_class op)
    // mpf_class sqrt (mpf_class op)
    // void mpf_class::swap (mpf_class& op)
    // void swap (mpf_class& op1, mpf_class& op2)
    // mpf_class trunc (mpf_class op)
    friend int sgn(const mpf_class &op);
    friend mpf_class sqrt(const mpf_class &op);
    friend mpf_class neg(const mpf_class &op);
    void swap(mpf_class &op) { mpf_swap(this->value, op.value); }
#if !defined ___GMPXX_STRICT_COMPATIBILITY___
    friend void swap(mpf_class &op1, mpf_class &op2) { mpf_swap(op1.value, op2.value); }
#endif
    friend mpf_class trunc(const mpf_class &op);

    // mp_bitcnt_t mpf_class::get_prec()
    // void mpf_class::set_prec (mp_bitcnt_t prec)
    // void mpf_class::set_prec_raw (mp_bitcnt_t prec)
    mp_bitcnt_t get_prec() const { return mpf_get_prec(value); }
    void set_prec(mp_bitcnt_t prec) { mpf_set_prec(value, prec); }
    void set_prec_raw(mp_bitcnt_t prec) { mpf_set_prec_raw(value, prec); }

    friend std::ostream &operator<<(std::ostream &os, const mpf_class &op);
    friend std::ostream &operator<<(std::ostream &os, const mpf_t op);

    friend std::istream &operator>>(std::istream &stream, mpf_class &op);
    friend std::istream &operator>>(std::istream &stream, mpf_t op);

    // gmpxx_mkII extension
    static mpf_class const_pi();
    static mpf_class const_e();
    static mpf_class const_log10();
    static mpf_class const_log2();
    static void reset_pi_cache();
    static void reset_e_cache();
    static void reset_log10_cache();
    static void reset_log2_cache();

    operator mpq_class() const;
    operator mpz_class() const;
    mpf_class &operator=(const mpz_class &other) {
        mpf_set_z(this->value, other.get_mpz_t());
        return *this;
    }
    mpf_class &operator=(const mpq_class &other) {
        mpf_set_q(this->value, other.get_mpq_t());
        return *this;
    }

    mpf_srcptr get_mpf_t() const { return value; }
    mpf_ptr get_mpf_t() { return value; }

  private:
    mpf_t value;
};
// casts
inline mpf_class::operator mpq_class() const { return mpq_class(this->get_mpf_t()); }
inline mpf_class::operator mpz_class() const { return mpz_class(this->get_mpf_t()); }
inline mpz_class::operator mpf_class() const { return mpf_class(this->get_mpz_t()); }
inline mpq_class::operator mpf_class() const { return mpf_class(this->get_mpq_t()); }
inline mpz_class::operator mpq_class() const { return mpq_class(this->get_mpz_t()); }
inline mpq_class::operator mpz_class() const { return mpz_class(this->get_mpq_t()); }

inline int preccmp(const mpf_class &op1, const mpf_class &op2) {
    mp_bitcnt_t prec1 = op1.get_prec();
    mp_bitcnt_t prec2 = op2.get_prec();
    if (prec1 > prec2) {
        return 1;
    } else if (prec1 == prec2) {
        return 0;
    } else {
        return -1;
    }
}
inline mpf_class &operator+=(mpf_class &lhs, const mpf_class &rhs) {
    mpf_add(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpf_class &operator-=(mpf_class &lhs, const mpf_class &rhs) {
    mpf_sub(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpf_class &operator*=(mpf_class &lhs, const mpf_class &rhs) {
    mpf_mul(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpf_class &operator/=(mpf_class &lhs, const mpf_class &rhs) {
    mpf_div(lhs.value, lhs.value, rhs.value);
    return lhs;
}
inline mpf_class operator+(const mpf_class &op1, const mpf_class &op2) {
#if !defined ___GMPXX_MKII_NOPRECCHANGE___
    mp_bitcnt_t max_prec = std::max(op1.get_prec(), op2.get_prec());
    mpf_class result(0, max_prec);
#else
    mpf_class result(0);
#endif
    mpf_add(result.get_mpf_t(), op1.get_mpf_t(), op2.get_mpf_t());
    return result;
}
inline mpf_class operator-(const mpf_class &op1, const mpf_class &op2) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class result(op1);
    mpf_sub(result.value, result.value, op2.value);
    return result;
#else
    int _preccmp = preccmp(op1, op2);
    if (_preccmp == 1) {
        mpf_class result(op1);
        mpf_sub(result.value, op1.value, op2.value);
        return result;
    } else {
        mpf_class result(op2);
        mpf_sub(result.value, op1.value, op2.value);
        return result;
    }
#endif
}
inline mpf_class operator*(const mpf_class &op1, const mpf_class &op2) {
#if !defined ___GMPXX_MKII_NOPRECCHANGE___
    mp_bitcnt_t max_prec = std::max(op1.get_prec(), op2.get_prec());
    mpf_class result(0, max_prec);
#else
    mpf_class result(0);
#endif
    mpf_mul(result.get_mpf_t(), op1.get_mpf_t(), op2.get_mpf_t());
    return result;
}
inline mpf_class operator/(const mpf_class &op1, const mpf_class &op2) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class result(op1);
    mpf_div(result.value, result.value, op2.value);
    return result;
#else
    int _preccmp = preccmp(op1, op2);
    if (_preccmp == 1) {
        mpf_class result(op1);
        mpf_div(result.value, op1.value, op2.value);
        return result;
    } else {
        mpf_class result(op2);
        mpf_div(result.value, op1.value, op2.value);
        return result;
    }
#endif
}
inline mpf_class operator+(const mpf_class &op) {
    mpf_class result(op.value);
    return result;
}
inline mpf_class operator-(const mpf_class &op) {
    mpf_class result(op.value);
    mpf_neg(result.value, op.value);
    return result;
}

// mpz_class cmp
inline int cmp(const mpz_class &op1, const mpz_class &op2) { return mpz_cmp(op1.get_mpz_t(), op2.get_mpz_t()); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const mpz_class &op1, T op2) { return mpz_cmp_ui(op1.get_mpz_t(), static_cast<unsigned long int>(op2)); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const T op1, mpz_class &op2) { return -mpz_cmp_ui(op2.get_mpz_t(), static_cast<unsigned long int>(op1)); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const mpz_class &op1, T op2) { return mpz_cmp_si(op1.get_mpz_t(), static_cast<signed long int>(op2)); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const T op1, mpz_class &op2) { return -mpz_cmp_si(op2.get_mpz_t(), static_cast<signed long int>(op1)); }
template <typename T> inline NON_INT_COND(T, int) cmp(const mpz_class &op1, T op2) { return mpz_cmp(op1.get_mpz_t(), mpz_class(op2).get_mpz_t()); }
template <typename T> inline NON_INT_COND(T, int) cmp(const T op1, mpz_class &op2) { return -mpz_cmp(op2.get_mpz_t(), mpz_class(op1).get_mpz_t()); }

// mpq_class cmp
inline int cmp(const mpq_class &op1, const mpq_class &op2) { return mpq_cmp(op1.get_mpq_t(), op2.get_mpq_t()); }
inline int cmp(const mpq_class &op1, const mpz_class &op2) { return mpq_cmp_z(op1.get_mpq_t(), op2.get_mpz_t()); }
inline int cmp(const mpz_class &op1, const mpq_class &op2) { return -mpq_cmp_z(op2.get_mpq_t(), op1.get_mpz_t()); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const mpq_class &op1, T op2) { return mpq_cmp(op1.get_mpq_t(), mpq_class(op2).get_mpq_t()); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const T op1, mpq_class &op2) { return -mpq_cmp(op2.get_mpq_t(), mpq_class(op1).get_mpq_t()); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const mpq_class &op1, T op2) { return mpq_cmp(op1.get_mpq_t(), mpq_class(op2).get_mpq_t()); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const T op1, mpq_class &op2) { return -mpq_cmp(op2.get_mpq_t(), mpq_class(op1).get_mpq_t()); }
template <typename T> inline NON_INT_COND(T, int) cmp(const mpq_class &op1, T op2) { return mpq_cmp(op1.get_mpq_t(), mpq_class(op2).get_mpq_t()); }
template <typename T> inline NON_INT_COND(T, int) cmp(const T op1, mpq_class &op2) { return -mpq_cmp(op2.get_mpq_t(), mpq_class(op1).get_mpq_t()); }

// mpf_class cmp
inline int cmp(const mpf_class &op1, const mpf_class &op2) { return mpf_cmp(op1.get_mpf_t(), op2.get_mpf_t()); }
inline int cmp(const mpf_class &op1, const mpq_class &op2) { return mpf_cmp(op1.get_mpf_t(), mpf_class(op2).get_mpf_t()); }
inline int cmp(const mpq_class &op1, const mpf_class &op2) { return mpf_cmp(mpf_class(op1).get_mpf_t(), op2.get_mpf_t()); }
inline int cmp(const mpf_class &op1, const mpz_class &op2) { return mpf_cmp_z(op1.get_mpf_t(), op2.get_mpz_t()); }
inline int cmp(const mpz_class &op1, const mpf_class &op2) { return -mpf_cmp_z(op2.get_mpf_t(), op1.get_mpz_t()); }
inline int cmp(const mpf_class &op1, const double op2) { return mpf_cmp_d(op1.get_mpf_t(), op2); }
inline int cmp(const double op1, const mpf_class &op2) { return -mpf_cmp_d(op2.get_mpf_t(), op1); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const mpf_class &op1, T op2) { return mpf_cmp_ui(op1.get_mpf_t(), static_cast<unsigned long int>(op2)); }
template <typename T> inline UNSIGNED_INT_COND(T, int) cmp(const T op1, mpf_class &op2) { return -mpf_cmp_ui(op2.get_mpf_t(), static_cast<unsigned long int>(op1)); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const mpf_class &op1, T op2) { return mpf_cmp_si(op1.get_mpf_t(), static_cast<signed long int>(op2)); }
template <typename T> inline SIGNED_INT_COND(T, int) cmp(const T op1, mpf_class &op2) { return -mpf_cmp_si(op2.get_mpf_t(), static_cast<signed long int>(op1)); }
template <typename T> inline NON_INT_COND(T, int) cmp(const mpf_class &op1, T op2) { return mpf_cmp(op1.get_mpf_t(), mpf_class(op2).get_mpf_t()); }
template <typename T> inline NON_INT_COND(T, int) cmp(const T op1, mpf_class &op2) { return -mpf_cmp(op2.get_mpf_t(), mpf_class(op1).get_mpf_t()); }

// implimentation of mpf_class operators
// improvements can be done using mpf_XXX_ui (note that they are not _si)
template <typename T> inline SIGNED_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs) {
    mpf_class _rhs(rhs);
    if (rhs >= 0) {
        mpf_add_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    } else {
        mpf_sub_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
    }
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    result += op2;
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2) { return op2 + op1; }
template <typename T> inline SIGNED_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs) {
    if (rhs >= 0) {
        mpf_sub_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    } else {
        mpf_add_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
    }
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    if (op2 >= 0) {
        mpf_sub_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    } else {
        mpf_add_ui(result.value, result.value, static_cast<unsigned long int>(-op2));
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2) {
    mpf_class result(op2);
    if (op1 >= 0) {
        mpf_ui_sub(result.value, static_cast<unsigned long int>(op1), op2.value);
    } else {
        mpf_add_ui(result.value, op2.value, static_cast<unsigned long int>(-op1));
        mpf_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs) {
    if (rhs >= 0) {
        mpf_mul_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    } else {
        mpf_mul_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
        mpf_neg(lhs.value, lhs.value);
    }
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    if (op2 >= 0) {
        mpf_mul_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    } else {
        mpf_mul_ui(result.value, result.value, static_cast<unsigned long int>(-op2));
        mpf_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2) { return op2 * op1; }
template <typename T> inline SIGNED_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs) {
    if (rhs >= 0) {
        mpf_div_ui(lhs.value, lhs.value, static_cast<unsigned long int>(rhs));
    } else {
        mpf_div_ui(lhs.value, lhs.value, static_cast<unsigned long int>(-rhs));
        mpf_neg(lhs.value, lhs.value);
    }
    return lhs;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    if (op2 >= 0) {
        mpf_div_ui(result.value, result.value, static_cast<unsigned long int>(op2));
    } else {
        mpf_div_ui(result.value, result.value, static_cast<unsigned long int>(-op2));
        mpf_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline SIGNED_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2) {
    mpf_class result(op2);
    if (op1 >= 0) {
        mpf_ui_div(result.value, static_cast<unsigned long int>(op1), result.value);
    } else {
        mpf_ui_div(result.value, static_cast<unsigned long int>(-op1), result.value);
        mpf_neg(result.value, result.value);
    }
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs) {
    mpf_add_ui(lhs.value, lhs.value, rhs);
    return lhs;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    mpf_add_ui(result.value, result.value, op2);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2) { return op2 + op1; }
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs) {
    mpf_sub_ui(lhs.value, lhs.value, rhs);
    return lhs;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    mpf_sub_ui(result.value, result.value, op2);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2) {
    mpf_class result(op2);
    mpf_ui_sub(result.value, op1, result.value);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs) {
    mpf_mul_ui(lhs.value, lhs.value, rhs);
    return lhs;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    mpf_mul_ui(result.value, result.value, op2);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2) { return op2 * op1; }
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs) {
    mpf_div_ui(lhs.value, lhs.value, rhs);
    return lhs;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    mpf_div_ui(result.value, op1.value, op2);
    return result;
}
template <typename T> inline UNSIGNED_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2) {
    mpf_class result(op2);
    mpf_ui_div(result.value, op1, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class &) operator+=(mpf_class &lhs, const T rhs) {
    mpf_class _rhs(rhs);
    mpf_add(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator+(const mpf_class &op1, const T op2) {
    mpf_class result(op1);
    result += op2;
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator+(const T op1, const mpf_class &op2) { return op2 + op1; }
template <typename T> inline NON_INT_COND(T, mpf_class &) operator-=(mpf_class &lhs, const T rhs) {
    mpf_class _rhs(rhs);
    mpf_sub(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator-(const mpf_class &op1, const T op2) {
    mpf_class result(op2);
    mpf_sub(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator-(const T op1, const mpf_class &op2) {
    mpf_class result(op1);
    mpf_sub(result.value, result.value, op2.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class &) operator*=(mpf_class &lhs, const T rhs) {
    mpf_class _rhs(rhs);
    mpf_mul(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator*(const mpf_class &op1, const T op2) {
    mpf_class result(op2);
    mpf_mul(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator*(const T op1, const mpf_class &op2) { return op2 * op1; }
template <typename T> inline NON_INT_COND(T, mpf_class &) operator/=(mpf_class &lhs, const T rhs) {
    mpf_class _rhs(rhs);
    mpf_div(lhs.value, lhs.value, _rhs.value);
    return lhs;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator/(const mpf_class &op1, const T op2) {
    mpf_class result(op2);
    mpf_div(result.value, op1.value, result.value);
    return result;
}
template <typename T> inline NON_INT_COND(T, mpf_class) operator/(const T op1, const mpf_class &op2) {
    mpf_class result(op1);
    mpf_div(result.value, result.value, op2.value);
    return result;
}
inline mpf_class &operator+=(mpf_class &lhs, const mpz_class &rhs) {
    mpf_class temp = mpf_class(rhs);
    lhs += temp;
    return lhs;
}
inline mpf_class &operator-=(mpf_class &lhs, const mpz_class &rhs) {
    mpf_class temp = mpf_class(rhs);
    lhs -= temp;
    return lhs;
}
inline mpf_class &operator*=(mpf_class &lhs, const mpz_class &rhs) {
    mpf_class temp = mpf_class(rhs);
    lhs *= temp;
    return lhs;
}
inline mpf_class &operator/=(mpf_class &lhs, const mpz_class &rhs) {
    mpf_class temp = mpf_class(rhs);
    lhs /= temp;
    return lhs;
}
inline mpf_class operator+(const mpf_class &op1, const mpz_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 += _op2;
    return _op1;
}
inline mpf_class operator+(const mpz_class &op1, const mpf_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op2 += _op1;
    return _op2;
}
inline mpf_class operator-(const mpf_class &op1, const mpz_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 -= _op2;
    return _op1;
}
inline mpf_class operator-(const mpz_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 -= _op2;
    return _op1;
}
inline mpf_class operator*(const mpf_class &op1, const mpz_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 *= _op2;
    return _op1;
}
inline mpf_class operator*(const mpz_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 *= _op2;
    return _op1;
}
inline mpf_class operator/(const mpf_class &op1, const mpz_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 /= _op2;
    return _op1;
}
inline mpf_class operator/(const mpz_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 /= _op2;
    return _op1;
}
inline mpf_class &operator+=(mpf_class &lhs, const mpq_class &rhs) {
    mpf_class temp(rhs);
    lhs += temp;
    return lhs;
}
inline mpf_class &operator-=(mpf_class &lhs, const mpq_class &rhs) {
    mpf_class temp(rhs);
    lhs -= temp;
    return lhs;
}

inline mpf_class &operator*=(mpf_class &lhs, const mpq_class &rhs) {
    mpf_class temp(rhs);
    lhs *= temp;
    return lhs;
}

inline mpf_class &operator/=(mpf_class &lhs, const mpq_class &rhs) {
    mpf_class temp(rhs);
    lhs /= temp;
    return lhs;
}
inline mpf_class operator+(const mpf_class &op1, const mpq_class &op2) {
    mpf_class result(op1);
    result += op2;
    return result;
}

inline mpf_class operator+(const mpq_class &op1, const mpf_class &op2) {
    mpf_class result(op2);
    result += op1;
    return result;
}

inline mpf_class operator-(const mpf_class &op1, const mpq_class &op2) {
    mpf_class result(op1);
    result -= op2;
    return result;
}
inline mpf_class operator-(const mpq_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 -= _op2;
    return _op1;
}
inline mpf_class operator*(const mpf_class &op1, const mpq_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 *= _op2;
    return _op1;
}
inline mpf_class operator*(const mpq_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 *= _op2;
    return _op1;
}
inline mpf_class operator/(const mpf_class &op1, const mpq_class &op2) {
    mpf_class _op1(op1);
    mpf_class _op2(op2);
    _op1 /= _op2;
    return _op1;
}
inline mpf_class operator/(const mpq_class &op1, const mpf_class &op2) {
    mpf_class _op1(op2); // 'op2' is used for correct precision initialization, furthur optimization may be possible for ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class _op2(op2);
    _op1 = op1;
    _op1 /= _op2;
    return _op1;
}

// elementary functions
inline mpf_class trunc(const mpf_class &op) {
    mpf_class rop(op);
    mpf_trunc(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class sqrt(const mpf_class &op) {
    mpf_class rop(op);
    mpf_sqrt(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class neg(const mpf_class &op) {
    mpf_class rop(op);
    mpf_neg(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class abs(const mpf_class &op) {
    mpf_class rop(op);
    mpf_abs(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class ceil(const mpf_class &op) {
    mpf_class rop(op);
    mpf_ceil(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class floor(const mpf_class &op) {
    mpf_class rop(op);
    mpf_floor(rop.value, op.get_mpf_t());
    return rop;
}
inline mpf_class hypot(const mpf_class &op1, const mpf_class &op2) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    mpf_class rop;
    rop = sqrt(op1 * op1 + op2 * op2);
    return rop;
#else
    int _preccmp = preccmp(op1, op2);
    if (_preccmp == 1) {
        mpf_class rop(op1);
        rop = sqrt(op1 * op1 + op2 * op2);
        return rop;
    } else {
        mpf_class rop(op2);
        rop = sqrt(op1 * op1 + op2 * op2);
        return rop;
    }
#endif
}
inline int sgn(const mpf_class &op) {
    int flag = mpf_sgn(op.get_mpf_t());
    return flag;
}
inline std::string mpf_to_base_string_default(const mpf_t value, int base, int flags, int width, int prec, char fill) {
    mp_exp_t exp;
    int effective_prec = (prec == 0) ? 6 : prec;
    char *base_cstr = mpf_get_str(nullptr, &exp, base, effective_prec, value);
    std::string base_str(base_cstr);
    free(base_cstr);

    bool is_showbase = flags & std::ios::showbase;
    bool is_showpoint = flags & std::ios::showpoint;
    bool is_uppercase = flags & std::ios::uppercase;
    std::string formatted_base;
    if (mpf_sgn(value) < 0) {
        base_str.erase(0, 1);
    }
    if (exp <= 0) {
        formatted_base = "0.";
        formatted_base.append(-exp, '0');
        formatted_base += base_str;
    } else if (size_t(exp) > base_str.length()) {
        formatted_base = base_str.substr(0, 1) + "." + base_str.substr(1);
        int adjusted_exp = exp - 1;
        std::string exp_str = adjusted_exp < base && adjusted_exp > -base ? "0" + std::to_string(adjusted_exp) : std::to_string(adjusted_exp);
        formatted_base += "e+" + exp_str; // can be minus?
    } else {
        formatted_base = base_str.substr(0, exp);
        if (exp < static_cast<mp_exp_t>(base_str.size())) {
            formatted_base += "." + base_str.substr(exp);
        }
    }
    if (is_showbase) {
        if (base == 16) {
            formatted_base.insert(0, "0x");
        } else if (base == 8 && mpf_sgn(value) != 0) {
            formatted_base.insert(0, "0");
        }
    }
    if (is_showpoint && formatted_base.find('.') == std::string::npos) {
        formatted_base += ".";
        while (formatted_base.length() < static_cast<size_t>(effective_prec + 1)) {
            formatted_base += '0';
        }
    }
    if (is_showpoint && base == 10 && formatted_base == "0.") {
        while (formatted_base.length() < static_cast<size_t>(effective_prec + 1)) {
            formatted_base += '0';
        }
    }
    if (mpf_sgn(value) < 0) {
        formatted_base.insert(0, "-");
    }
    if (width > static_cast<int>(formatted_base.size())) {
        std::streamsize padding_length = width - formatted_base.size();
        if (flags & std::ios_base::left) {
            formatted_base.append(padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 16 && formatted_base[0] == '0' && formatted_base[1] == 'x') { // Insert padding after the "0x"
            formatted_base.insert(2, padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 16 && formatted_base[0] == '-' && formatted_base[1] == '0' && formatted_base[2] == 'x') { // Insert padding after the "-0x"
            formatted_base.insert(3, padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 10) {
            size_t pos = 0;
            if (formatted_base[0] == '-' || formatted_base[0] == '+') {
                pos = 1;
            }
            formatted_base.insert(pos, padding_length, fill);
        } else {
            formatted_base.insert(0, padding_length, fill);
        }
    }
    if (!is_showpoint) {
        if (formatted_base.back() == '.') {
            formatted_base.erase(formatted_base.size() - 1);
        }
    }
    if (is_uppercase) {
        std::transform(formatted_base.begin(), formatted_base.end(), formatted_base.begin(), [](unsigned char c) { return std::toupper(c); });
    }
    return formatted_base;
}
inline mp_exp_t integraldigits_in_base(const mpf_t value, int base) {
    mp_exp_t exp;
    mpf_get_d_2exp(&exp, value);
    if (mpf_sgn(value) == 0) {
        return 1;
    }
    if (exp <= 0) {
        return 1;
    }
    int digits = 1;
    switch (base) {
    case 16:
        digits = (exp + 3) / 4;
        break;
    case 8:
        digits = (exp + 2) / 3;
        break;
    case 10:
        digits = static_cast<int>(std::floor(exp * 1.0 / std::log2(10)) + 1); // the number of digits can sometimes be larger than the actual digits
        break;
    case 2:
        digits = exp;
        break;
    default:
        return 0;
    }
    return digits;
}
inline std::string mpf_to_base_string_fixed(const mpf_t value, int base, int flags, int width, int prec, char fill) {
    mp_exp_t exp;
    int effective_prec = (prec == 0) ? 6 : prec;
    mp_exp_t digits = integraldigits_in_base(value, base);
    char *base_cstr = mpf_get_str(nullptr, &exp, base, digits + effective_prec, value);
    std::string base_str(base_cstr);
    free(base_cstr);
    bool is_showbase = flags & std::ios::showbase;
    bool is_showpoint = flags & std::ios::showpoint;
    bool is_uppercase = flags & std::ios::uppercase;
    std::string formatted_base;

    if (mpf_sgn(value) == 0) {
        if (base == 16) {
            formatted_base.insert(0, "0x0");
        } else if (base == 10) {
            formatted_base.insert(0, "0");
            if (prec == 0) {
                if (is_showpoint) {
                    formatted_base += ".";
                }
            } else {
                formatted_base += ".";
                formatted_base += std::string(effective_prec, '0');
            }
        } else {
            formatted_base.insert(0, "0");
        }
    } else {
        if (exp > 0) {
            formatted_base += base_str.substr(0, exp); // integer part
            if (effective_prec > 0 || is_showpoint) {
                formatted_base += ".";
                formatted_base += base_str.substr(exp, effective_prec); // fraction part
            }
        } else {
            formatted_base += "0.";
            std::string tobe_added = std::string(-exp, '0') + base_str;
            formatted_base += tobe_added.substr(0, effective_prec); // 0.000XXXX
        }
        if (is_showbase) {
            if (base == 16) {
                formatted_base.insert(0, "0x");
            } else if (base == 8 && formatted_base != "0.0") {
                formatted_base.insert(0, "0");
            }
        }
        if (prec != 0 && (!formatted_base.empty() && formatted_base.back() == '.')) {
            formatted_base += std::string(effective_prec, '0');
        }
    }
    if (width > static_cast<int>(formatted_base.size())) {
        std::streamsize padding_length = width - formatted_base.size();
        if (flags & std::ios_base::left) {
            formatted_base.append(padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 16 && formatted_base[0] == '0' && formatted_base[1] == 'x') { // Insert padding after the "0x"
            formatted_base.insert(2, padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 16 && formatted_base[0] == '-' && formatted_base[1] == '0' && formatted_base[2] == 'x') { // Insert padding after the "-0x"
            formatted_base.insert(3, padding_length, fill);
        } else if (flags & std::ios_base::internal && base == 10) {
            size_t pos = 0;
            if (formatted_base[0] == '-' || formatted_base[0] == '+') {
                pos = 1;
            }
            formatted_base.insert(pos, padding_length, fill);
        } else {
            formatted_base.insert(0, padding_length, fill);
        }
    }
    if (!is_showpoint) {
        if (formatted_base.back() == '.') {
            formatted_base.erase(formatted_base.size() - 1);
        }
    }
    if (is_uppercase) {
        std::transform(formatted_base.begin(), formatted_base.end(), formatted_base.begin(), [](unsigned char c) { return std::toupper(c); });
    }
    return formatted_base;
}
inline std::string mpf_to_base_string_scientific(const mpf_t value, int base, int flags, int width, int prec, char fill) {
    // TODO obtain correct # of digits in the given base, check rounding of the negative exp part
    mp_exp_t exp;
    int effective_prec = (prec == 0) ? 6 : prec;
    char *base_cstr = mpf_get_str(nullptr, &exp, base, effective_prec + 1, value);
    std::string base_str(base_cstr);
    free(base_cstr);
    bool is_showbase = flags & std::ios::showbase;
    bool is_uppercase = flags & std::ios::uppercase;
    std::string formatted_base;
    if (mpf_sgn(value) == 0) {
        formatted_base = "0." + std::string(effective_prec, '0');
        exp = 1;
    } else {
        if (base_str[0] == '-') {
            base_str.erase(0, 1);
            formatted_base += "-";
        }
        formatted_base += base_str[0];
        if (base_str.length() > 1) {
            formatted_base += "." + base_str.substr(1, effective_prec);
        } else {
            formatted_base += "." + std::string(effective_prec, '0');
        }
    }
    int _pad = 2;
    if (mpf_sgn(value) < 0)
        _pad++;
    if (formatted_base.length() < static_cast<size_t>(effective_prec + _pad)) {
        size_t padding_length = effective_prec + _pad - formatted_base.length();
        formatted_base.append(padding_length, '0');
    }
    int adjusted_exp = exp - 1;
    if (base == 16) {
        formatted_base += "@";
    } else if (base == 8 || base == 10) {
        formatted_base += "e";
    }
    if (adjusted_exp >= 0) {
        formatted_base += "+";
        formatted_base += (adjusted_exp < 10 ? "0" : "") + std::to_string(adjusted_exp);
    } else {
        formatted_base += "-";
        formatted_base += (-adjusted_exp < 10 ? "0" : "") + std::to_string(-adjusted_exp);
    }
    if (is_showbase) {
        if (base == 16) {
            if (formatted_base[0] == '-' || formatted_base[0] == '+') {
                formatted_base.insert(1, "0x");
            } else {
                formatted_base.insert(0, "0x");
            }
        } else if (base == 8) {
            if (formatted_base[0] == '-' || formatted_base[0] == '+') {
                formatted_base.insert(1, "0");
            } else {
                formatted_base.insert(0, "0");
            }
        }
    }
    if (static_cast<int>(formatted_base.size()) < width) {
        int padding_length = width - formatted_base.size();
        if (flags & std::ios_base::left) {
            formatted_base.append(padding_length, fill);
        } else if (flags & std::ios_base::internal && formatted_base[0] == '-') {
            size_t pos = formatted_base.find_first_not_of('-');
            formatted_base.insert(pos, padding_length, fill);
        } else {
            formatted_base.insert(0, padding_length, fill);
        }
    }
    if (is_uppercase) {
        std::transform(formatted_base.begin(), formatted_base.end(), formatted_base.begin(), [](unsigned char c) { return std::toupper(c); });
    }
    return formatted_base;
}
inline void print_mpf(std::ostream &os, const mpf_t op) {
    std::ios_base::fmtflags flags = os.flags();
    std::streamsize prec = os.precision();
    std::streamsize width = os.width();
    bool is_hex = flags & std::ios::hex;
    bool is_oct = flags & std::ios::oct;
    bool is_dec = flags & std::ios::dec;
    bool is_fixed = flags & std::ios::fixed;
    bool is_scientific = flags & std::ios::scientific;
    char fill = os.fill();
    char *str = nullptr;

    int base = 10;
    if (is_hex) {
        base = 16;
    } else if (is_dec) {
        base = 10;
    } else if (is_oct) {
        base = 8;
    }
    std::string format;
    std::string base_string;
    if (is_fixed) {
        base_string = mpf_to_base_string_fixed(op, base, flags, width, prec, fill);
    } else if (is_scientific) {
        base_string = mpf_to_base_string_scientific(op, base, flags, width, prec, fill);
    } else {
        base_string = mpf_to_base_string_default(op, base, flags, width, prec, fill);
    }
    str = strdup(base_string.c_str());
    std::string s(str);
    free(str);
    if (flags & std::ios::showpos && mpf_sgn(op) >= 0) {
        s.insert(0, "+");
    }
    std::streamsize len = s.length();
    if (len < width) {
        std::streamsize padding_length = width - len;
        if (flags & std::ios::left) {
            s.append(padding_length, fill);
        } else if (flags & std::ios::internal && s[0] == '-') {
            size_t pos = s.find_first_not_of('-');
            s.insert(pos, padding_length, fill);
        } else {
            s.insert(0, padding_length, fill);
        }
    }
    os << s;
    os.width(0);
}
inline std::ostream &operator<<(std::ostream &os, const mpf_class &op) {
    print_mpf(os, op.get_mpf_t());
    return os;
}
inline std::ostream &operator<<(std::ostream &os, const mpf_t &op) {
    print_mpf(os, op);
    return os;
}
inline bool is_valid_number_char(char ch) { return std::isxdigit(ch) || ch == '.' || ch == 'p' || ch == 'P' || ch == '-' || ch == '+'; }
inline void print_format_flags(std::ios_base::fmtflags flags) {
    std::cout << "Current Format Flags:" << std::endl;
    if (flags & std::ios_base::dec)
        std::cout << "dec ";
    if (flags & std::ios_base::oct)
        std::cout << "oct ";
    if (flags & std::ios_base::hex)
        std::cout << "hex ";
    if (flags & std::ios_base::scientific)
        std::cout << "scientific ";
    if (flags & std::ios_base::fixed)
        std::cout << "fixed ";
    if (flags & std::ios_base::boolalpha)
        std::cout << "boolalpha ";
    if (flags & std::ios_base::showbase)
        std::cout << "showbase ";
    if (flags & std::ios_base::showpoint)
        std::cout << "showpoint ";
    if (flags & std::ios_base::showpos)
        std::cout << "showpos ";
    if (flags & std::ios_base::skipws)
        std::cout << "skipws ";
    if (flags & std::ios_base::unitbuf)
        std::cout << "unitbuf ";
    if (flags & std::ios_base::uppercase)
        std::cout << "uppercase ";
    if (flags & std::ios_base::adjustfield)
        std::cout << "adjustfield ";
    if (flags & std::ios_base::basefield)
        std::cout << "basefield ";
    if (flags & std::ios_base::floatfield)
        std::cout << "floatfield ";
    std::cout << std::endl;
}
inline std::istream &read_mpf_from_stream(std::istream &stream, mpf_t op) {
    std::ios_base::fmtflags current_flags = stream.flags();
    if (current_flags & std::ios_base::oct || current_flags & std::ios_base::hex) {
        throw std::runtime_error("Unsupported number base for mpf_t");
    }
    char ch;
    std::string number;
    bool negative = false;
    bool is_space = false;
    int base = 10;
    int counter = 0;
    while (stream >> ch && isspace(ch)) {
        is_space = true;
        counter++;
    }
    if (!(current_flags & std::ios::skipws) && is_space == true) {
        for (int i = 0; i <= counter; i++)
            stream.unget();
        stream.setstate(std::ios::failbit);
        return stream;
    }
    if (ch == '+' || ch == '-') {
        negative = (ch == '-');
        if (!stream.get(ch)) {
            stream.setstate(std::ios::failbit);
            return stream;
        }
        if (ch == '+' || ch == '-') {
            stream.unget();
            stream.setstate(std::ios::failbit);
            return stream;
        }
    }
    if (ch == '.') {
        if (!stream.get(ch)) {
            stream.setstate(std::ios::failbit);
            return stream;
        }
        if (ch == 'e' || ch == 'E') {
            stream.unget();
            stream.setstate(std::ios::failbit);
            return stream;
        }
    }
    if (ch == 'e' || ch == 'E') {
        stream.unget();
        stream.setstate(std::ios::failbit);
        return stream;
    }
    if (!std::isdigit(ch) && ch != '.') {
        stream.setstate(std::ios::failbit);
        return stream;
    }
    while (is_valid_number_char(ch)) {
        number += ch;
        if (!stream.get(ch))
            break;
    }
    // The following code detects invalid +-, -+, --, ++
    std::size_t pos_plus_minus = number.find("+-");
    std::size_t pos_minus_plus = number.find("-+");
    std::size_t pos_minus_minus = number.find("--");
    std::size_t pos_plus_plus = number.find("++");
    std::size_t invalid_pos = std::string::npos;
    if (pos_plus_minus != std::string::npos) {
        invalid_pos = pos_plus_minus;
    } else if (pos_minus_plus != std::string::npos) {
        invalid_pos = pos_minus_plus;
    } else if (pos_minus_minus != std::string::npos) {
        invalid_pos = pos_minus_minus;
    } else if (pos_plus_plus != std::string::npos) {
        invalid_pos = pos_plus_plus;
    }
    if (invalid_pos != std::string::npos) {
        if (stream.eof())
            stream.clear();
        for (std::size_t i = number.size(); i > invalid_pos + 1; --i) {
            stream.unget();
        }
        stream.setstate(std::ios::failbit);
        return stream;
    }
    int ret = mpf_set_str(op, number.c_str(), base);
    if (ret != 0) {
        stream.setstate(stream.rdstate() & ~std::ios::goodbit);
        stream.setstate(std::ios::failbit);
    } else {
        stream.clear(stream.rdstate() & ~std::ios::failbit);
        stream.setstate(std::ios::goodbit);
    }
    if (negative) {
        mpf_neg(op, op);
    }
    return stream;
}
inline std::istream &operator>>(std::istream &stream, mpf_t op) { return read_mpf_from_stream(stream, op); }
inline std::istream &operator>>(std::istream &stream, mpf_class &op) { return read_mpf_from_stream(stream, op.get_mpf_t()); }
inline mpf_class const_pi() {
    static bool calculated = false;
    static mp_bitcnt_t calculated_pi_precision = 0;
    mp_bitcnt_t _default_prec = mpf_get_default_prec();
    if (!calculated || (calculated && calculated_pi_precision != _default_prec)) {
        calculated_pi_precision = mpf_get_default_prec();
        // calculating approximate pi using arithmetic-geometric mean
        mpf_class zero(0.0, _default_prec);
        mpf_class quarter(0.25, _default_prec);
        mpf_class one(1.0, _default_prec);
        mpf_class two(2.0, _default_prec);
        mpf_class four(4.0, _default_prec);
        mpf_class a(one, _default_prec), b(one / sqrt(two), _default_prec), t(quarter, _default_prec), p(one, _default_prec);
        mpf_class a_next(zero, _default_prec), b_next(zero, _default_prec), t_next(zero, _default_prec), tmp_pi(zero, _default_prec), pi_previous(zero, _default_prec);
        caches<>::pi_cached.set_prec(_default_prec);
        caches<>::pi_cached = zero;

        bool converged = false;
        int iteration = 0;

        mpf_class epsilon = one;
        epsilon.div_2exp(_default_prec);

        while (!converged) {
            iteration++;
            a_next = (a + b) / two;
            b_next = sqrt(a * b);
            t_next = t - p * (a - a_next) * (a - a_next);
            p = two * p;

            // Update values for the next iteration
            a = a_next;
            b = b_next;
            t = t_next;

            // Calculate pi
            pi_previous = tmp_pi;
            tmp_pi = (a + b) * (a + b) / (four * t);
            // Check for convergence
            if (abs(tmp_pi - pi_previous) < epsilon) {
                converged = true;
            }
        }
        calculated = true;
        caches<>::pi_cached = tmp_pi;
    } else {
        //      std::cout << "pi cached\n";
    }
    return caches<>::pi_cached;
}
inline mpf_class const_pi(mp_bitcnt_t req_precision) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // calculating approximate pi using arithmetic-geometric mean
    mpf_class zero(0.0, req_precision);
    mpf_class quarter(0.25, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class four(4.0, req_precision);

    mpf_class calculated_pi(zero, req_precision);
    mpf_class a(one, req_precision), b(one / sqrt(two), req_precision), t(quarter, req_precision), p(one, req_precision);
    mpf_class a_next(zero, req_precision), b_next(zero, req_precision), t_next(zero, req_precision), tmp_pi(zero, req_precision), pi_previous(zero, req_precision);
    mpf_class epsilon(zero, req_precision), tmp(zero, req_precision);

    assert(calculated_pi.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(t.get_prec() == req_precision);
    assert(p.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(t_next.get_prec() == req_precision);
    assert(tmp_pi.get_prec() == req_precision);
    assert(pi_previous.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);

    bool converged = false;
    int iteration = 0;

    epsilon = one;
    epsilon.div_2exp(req_precision);
    while (!converged) {
        iteration++;
        a_next = (a + b) / two;
        b_next = sqrt(a * b);
        t_next = t - p * (a - a_next) * (a - a_next);
        p = two * p;

        // Update values for the next iteration
        a = a_next;
        b = b_next;
        t = t_next;

        // Calculate pi
        pi_previous = tmp_pi;
        tmp_pi = (a + b) * (a + b) / (four * t);

        // Check for convergence
        tmp = abs(tmp_pi - pi_previous);
        if (tmp < epsilon) {
            converged = true;
        }
    }
    calculated_pi = tmp_pi;
    assert(calculated_pi.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(t.get_prec() == req_precision);
    assert(p.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(t_next.get_prec() == req_precision);
    assert(tmp_pi.get_prec() == req_precision);
    assert(pi_previous.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);

    return calculated_pi;
}

inline mpf_class const_log2() {
    static mpf_class log2_cached;
    static bool calculated = false;
    static mp_bitcnt_t calculated_log2_precision = 0;
    mp_bitcnt_t _default_prec = mpf_get_default_prec();

    if (!calculated || (calculated && calculated_log2_precision != _default_prec)) {
        log2_cached = mpf_class();
        calculated_log2_precision = mpf_get_default_prec();
        // calculating approximate log2 using arithmetic-geometric mean
        mpf_class zero(0.0);
        mpf_class one(1.0);
        mpf_class two(2.0);
        mpf_class a(one);
        mpf_class epsilon = one;
        epsilon.div_2exp((_default_prec / 2) - 2);

        log2_cached.set_prec(_default_prec);
        log2_cached = zero;

        mpf_class b = epsilon;
        mpf_class sum = one;
        mpf_class a_next, b_next;
        mpf_class tmp;

        bool converged = false;

        while (!converged) {
            a_next = (a + b) / two;
            b_next = sqrt(a * b);

            // Check for convergence
            if (abs(a - b) < epsilon) {
                converged = true;
            }
            a = a_next;
            b = b_next;
        }
        log2_cached = const_pi() / (mpf_class(_default_prec) * a);
        calculated = true;
    }
    return log2_cached;
}
inline mpf_class const_log2(mp_bitcnt_t req_precision) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);

    mpf_class log2(zero);
    mpf_class epsilon(one), tmp(zero);
    mpf_class a(one), b(one);
    mpf_class a_next(zero), b_next(zero);
    mpf_class sum(one);

    bool converged = false;

    // calculating approximate log2 using arithmetic-geometric mean
    b.div_2exp((req_precision / 2) - 2);
    epsilon.div_2exp(req_precision);

    assert(log2.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);
    assert(sum.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(one.get_prec() == req_precision);
    assert(two.get_prec() == req_precision);

    while (!converged) {
        a_next = (a + b) / two;
        b_next = sqrt(a * b);
        assert(b_next.get_prec() == req_precision);

        // Check for convergence
        if (abs(a - b) < epsilon) {
            converged = true;
        }
        a = a_next;
        b = b_next;
    }
    log2 = const_pi(req_precision) / (mpf_class(req_precision, req_precision) * a);

    assert(const_pi(req_precision).get_prec() == req_precision);
    assert(mpf_class(req_precision, req_precision).get_prec() == req_precision);
    assert(log2.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);
    assert(sum.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(one.get_prec() == req_precision);
    assert(two.get_prec() == req_precision);

    return log2;
}
inline mpf_class log(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class four(4.0, req_precision);

    mpf_class _log(zero);
    mpf_class epsilon(one), tmp(zero);
    mpf_class a(one), b(one);
    mpf_class a_next(zero), b_next(zero);
    mpf_class s(one);
    mpf_class _pi(const_pi(req_precision));
    mpf_class _log2(const_log2(req_precision));
    mp_exp_t m;
    bool converged = false;

    assert(_log.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(s.get_prec() == req_precision);
    assert(_pi.get_prec() == req_precision);
    assert(_log2.get_prec() == req_precision);

    // calculating approximate log2 using arithmetic-geometric mean
    b = one;
    b.mul_2exp(req_precision / 2);
    s = b / x;
    mpf_get_d_2exp(&m, s.get_mpf_t());

    b = one;
    b.mul_2exp(m);
    s = x * b;

    b = four / s;
    epsilon.div_2exp(req_precision);
    int counter = 0;
    while (!converged) {
        counter++;
        a_next = (a + b) / two;
        b_next = sqrt(a * b);

        // Check for convergence
        if (abs(a - b) < epsilon) {
            converged = true;
        }
        a = a_next;
        b = b_next;
    }
    _log = _pi / (two * b) - m * _log2;

    assert(_log.get_prec() == req_precision);
    assert(epsilon.get_prec() == req_precision);
    assert(tmp.get_prec() == req_precision);
    assert(a.get_prec() == req_precision);
    assert(b.get_prec() == req_precision);
    assert(a_next.get_prec() == req_precision);
    assert(b_next.get_prec() == req_precision);
    assert(s.get_prec() == req_precision);
    assert(_pi.get_prec() == req_precision);
    assert(_log2.get_prec() == req_precision);

    return _log;
}
inline mpf_class exp(const mpf_class &x) {
    // https://www.mpfr.org/algorithms.pdf section 4.4
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class _exp(one);
    mpf_class _x(x);
    mpf_class r(zero);
    mpf_class _pi(const_pi(req_precision));
    mpf_class _log2(const_log2(req_precision));
    mp_exp_t k = 0, l = 0, n = 0;

    if (x < zero)
        _x = -_x; // avoid cancellation of significant digits
    // calculating approximate exp
    // taking modulo of log2
    mpf_get_d_2exp(&k, _x.get_mpf_t());
    if (k > 0) {
        _x.div_2exp(k);    // 0.5 <= |x| < 1
        _log2.div_2exp(k); // log2/2 = 0.346574
        n = floor(_x / _log2).get_si();
        r = _x - n * _log2;
        l = req_precision / k;
    } else {
        k = 0;
        l = req_precision;
        r = _x;
        n = 0;
    }
    for (int i = l; i > 0; i--) {
        _exp = one + ((r * _exp) / mpf_class(i, req_precision));
    }
    for (int i = 0; i < k; i++) {
        _exp = _exp * _exp;
    }
    if (n > 0)
        _exp.mul_2exp(n);
    if (n < 0)
        _exp.div_2exp(-n);
    if (x < zero)
        _exp = one / _exp; // avoid cancellation of significant digits
    return _exp;
}
inline mpf_class mpf_remainder(const mpf_class &x, const mpf_class &y, mpz_class *quotient_out = nullptr) {
    mpf_class quotient = x / y;
    mpz_class int_quotient(quotient);
    if (quotient_out) {
        *quotient_out = int_quotient;
    }
    mpf_class remainder = x - int_quotient * y;
    if (remainder < 0) {
        remainder += y;
        if (quotient_out) {
            (*quotient_out)--;
        }
    }
    return remainder;
}
// Naive Taylor expansion version. It generates a very long series.
inline mpf_class cos_taylor_naive(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // Constants and variables
    mpf_class _PI(0.0, req_precision);
    mpf_class two_pi(0.0, req_precision);
    mpf_class pi_over_4(0.0, req_precision);
    mpf_class term(0.0, req_precision);
    mpf_class cosx(0.0, req_precision);
    mpf_class cosx_prev(0.0, req_precision);
    mpf_class x_squared(0.0, req_precision);
    mpf_class x_reduced(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class four(4.0, req_precision);
    mpf_class n(0.0, req_precision);
    mpf_class epsilon(0.0, req_precision);
    int sign = -1;
    int symm_sign = 1;
    // Setting some constants
    _PI = const_pi(req_precision);
    two_pi = two * _PI;
    pi_over_4 = _PI / four;
    epsilon.set_epsilon();
    // cos(-x) = cos(x)
    x_reduced = x;
    if (x_reduced < 0) {
        x_reduced = -x_reduced;
    }
    // Reduce x to [-pi, pi)
    x_reduced = mpf_remainder(x_reduced + _PI, two_pi);
    if (x_reduced < 0)
        x_reduced += two_pi;
    x_reduced -= _PI;
    // Furthur reduce x to  [-pi/2, pi/2)
    if (x_reduced > _PI / 2) {
        x_reduced = _PI - x;
        symm_sign = -1;
    } else if (x_reduced < -_PI / 2) {
        x_reduced = -_PI - x_reduced;
        symm_sign = -1;
    }
    // Calculate cos(x) using Taylor series
    x_squared = x_reduced * x_reduced;
    term = one;
    cosx = one;
    cosx_prev = one;
    for (mp_bitcnt_t _n = 1; _n < req_precision; _n++) {
        n = mpf_class(_n, req_precision);
        term *= x_reduced * x_reduced / (2 * n * (2 * n - 1));
        cosx += sign * term;
        if (abs(cosx_prev - cosx) < epsilon)
            break;
        cosx_prev = cosx;
        sign *= -1;
    }
    return cosx * symm_sign;
}
// Internal use only. calculate cos(x) where x is [-pi/2, pi/2).
inline mpf_class cos_taylor_reduced(const mpf_class &x, bool addprec = false) {
    mp_bitcnt_t _req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(_req_precision == mpf_get_default_prec());
#endif
    // Calculate cos(x) using Taylor series
    mp_bitcnt_t k = std::floor(std::sqrt(_req_precision / 2));
    // We need some additional precision for successive applications of the double-angle formula
    mp_bitcnt_t additional_precision = 0;
    if (k % 64 != 0) {
        additional_precision = ((k / 64) + 1) * 64;
    }
    mp_bitcnt_t req_precision = _req_precision + additional_precision;
    // Constants and variables
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class r(0.0, req_precision);
    mpf_class s(0.0, req_precision);
    mpf_class t(0.0, req_precision);
    mpf_class l(0.0, req_precision);
    mpf_class epsilon(1.0, req_precision);
    mpf_class _s(0.0, _req_precision);
    r = x;
    r = r * r;
    r.div_2exp((k * 2));
    epsilon.div_2exp(req_precision);
    s = zero;
    t = one;
    l = one;
    int counter = 0;
    for (mp_bitcnt_t _l = 1; _l < _req_precision; _l++) {
        t *= r;
        t /= ((two * l - one) * (two * l));
        if (t < epsilon) {
            break;
        }
        if (_l % 2 == 1) {
            s -= t;
        } else {
            s += t;
        }
        l = l + 1;
        counter++;
    }
    s += one;
    for (mp_bitcnt_t i = 0; i < k; i++) {
        s *= two * s;
        s -= one;
    }
    if (addprec == true)
        return s;
    _s = s; // reduce the precision, may not be necessary but to be sure.
    return _s;
}
inline mpf_class cos_taylor(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class pi(0.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class two_pi(0.0, req_precision);
    mpf_class pi_over_2(0.0, req_precision);
    mpf_class x_reduced(0.0, req_precision);
    mpf_class cosx(0.0, req_precision);
    int symm_sign = 1;
    // Setting some constants
    pi = const_pi(req_precision);
    two_pi = two * pi;
    pi_over_2 = pi / two;
    // cos(-x) = cos(x)
    x_reduced = x;
    if (x_reduced < 0) {
        x_reduced = -x_reduced;
    }
    // Reduce x to [-pi, pi)
    x_reduced = mpf_remainder(x_reduced + pi, two_pi);
    if (x_reduced < 0)
        x_reduced += two_pi;
    x_reduced -= pi;
    // Furthur reduce x to  [-pi/2, pi/2)
    if (x_reduced > pi_over_2) {
        x_reduced = pi - x;
        symm_sign = -1;
    } else if (x_reduced < -pi_over_2) {
        x_reduced = -pi - x_reduced;
        symm_sign = -1;
    }
    cosx = cos_taylor_reduced(x_reduced) * symm_sign;
    return cosx;
}
inline mpf_class cos(const mpf_class &x) { return cos_taylor(x); }
// mpf_class cos(const mpf_class &x) { return cos_taylor_naive(x); }
//  Naive Taylor expansion version. It generates a very long series.
inline mpf_class sin_taylor_naive(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // Constants and variables
    mpf_class _PI(0.0, req_precision);
    mpf_class two_pi(0.0, req_precision);
    mpf_class pi_over_2(0.0, req_precision);
    mpf_class term(0.0, req_precision);
    mpf_class sinx(0.0, req_precision);
    mpf_class sinx_prev(0.0, req_precision);
    mpf_class x_squared(0.0, req_precision);
    mpf_class x_reduced(0.0, req_precision);
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class three(3.0, req_precision);
    mpf_class four(4.0, req_precision);
    mpf_class n(0.0, req_precision);
    mpf_class epsilon(0.0, req_precision);
    int symm_sign = 1;
    // Setting some constants
    _PI = const_pi(req_precision);
    two_pi = two * _PI;
    pi_over_2 = _PI / two;
    epsilon.set_epsilon();
    // sin(-x) = -sin(x)
    x_reduced = x;
    if (x_reduced < 0) {
        x_reduced = -x_reduced;
        symm_sign = -1;
    }
    // Reduce x to [0, 2pi)
    x_reduced = mpf_remainder(x_reduced, two_pi);
    // Furthur reduce x to [0, pi/2)
    if ((pi_over_2 < x_reduced) && (x_reduced <= _PI)) {
        x_reduced = _PI - x_reduced;
    }
    if ((_PI < x_reduced) && (x_reduced <= three * two_pi)) {
        x_reduced = three * two_pi - x_reduced;
        symm_sign *= -1;
    }
    if ((three * two_pi < x_reduced) && (x_reduced <= two_pi)) {
        x_reduced = two_pi - x_reduced;
        symm_sign *= -1;
    }
    // Calculate sin(x) using Taylor series
    term = x_reduced;
    sinx = zero;
    sinx_prev = zero;
    for (mp_bitcnt_t _n = 1; _n < req_precision; _n++) {
        n = mpf_class(_n, req_precision);
        sinx += term;
        term *= -x_reduced * x_reduced / (two * n * (two * n + one));
        if (abs(sinx_prev - sinx) < epsilon)
            break;
        sinx_prev = sinx;
    }
    return sinx * symm_sign;
}
// assume x is in the interval [0, pi/2)
inline mpf_class sinx_from_cos_internal(const mpf_class &x, bool addprec = false) {
    mp_bitcnt_t _req_precision = x.get_prec();
    mpf_class _s(0.0, _req_precision);
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(_req_precision == mpf_get_default_prec());
#endif
    mp_bitcnt_t k = std::floor(std::sqrt(_req_precision / 2));
    // We need some additional precision for successive applications of the double-angle formula
    mp_bitcnt_t additional_precision = 0;
    if (k % 64 != 0) {
        additional_precision = ((k / 64) + 1) * 64;
    }
    mp_bitcnt_t req_precision = _req_precision + additional_precision;
    mpf_class c(0.0, req_precision);
    mpf_class t(0.0, req_precision);
    mpf_class s(0.0, req_precision);
    mpf_class u(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    c = cos_taylor_reduced(x, true);
    t = c * c;
    u = one - t;
    s = sqrt(u);
    if (addprec)
        return s;
    else
        _s = s;
    return _s;
}
inline mpf_class sin_from_cos(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // Constants and variables
    mpf_class _PI(0.0, req_precision);
    mpf_class two_pi(0.0, req_precision);
    mpf_class pi_over_2(0.0, req_precision);
    mpf_class x_reduced(0.0, req_precision);
    mpf_class sinx(0.0, req_precision);
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class three(3.0, req_precision);
    mpf_class n(0.0, req_precision);
    int symm_sign = 1;
    // Setting some constants
    _PI = const_pi(req_precision);
    two_pi = two * _PI;
    pi_over_2 = _PI / two;
    // sin(-x) = -sin(x)
    x_reduced = x;
    if (x_reduced < 0) {
        x_reduced = -x_reduced;
        symm_sign = -1;
    }
    // Reduce x to [0, 2pi)
    x_reduced = mpf_remainder(x_reduced, two_pi);
    // Furthur reduce x to [0, pi/2)
    if ((pi_over_2 < x_reduced) && (x_reduced <= _PI)) {
        x_reduced = _PI - x_reduced;
    }
    if ((_PI < x_reduced) && (x_reduced <= three * two_pi)) {
        x_reduced = three * two_pi - x_reduced;
        symm_sign *= -1;
    }
    if ((three * two_pi < x_reduced) && (x_reduced <= two_pi)) {
        x_reduced = two_pi - x_reduced;
        symm_sign *= -1;
    }
    sinx = sinx_from_cos_internal(x_reduced);
    return sinx * symm_sign;
}
inline mpf_class sin(const mpf_class &x) { return sin_from_cos(x); }
inline mpf_class tan_from_sin_cos(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // Constants and variables
    mpf_class _PI(0.0, req_precision);
    mpf_class two_pi(0.0, req_precision);
    mpf_class pi_over_2(0.0, req_precision);
    mpf_class x_reduced(0.0, req_precision);
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class three(3.0, req_precision);
    mpf_class n(0.0, req_precision);
    int symm_sign = 1;
    // Setting some constants
    _PI = const_pi(req_precision);
    two_pi = two * _PI;
    pi_over_2 = _PI / two;
    // tan(-x) = -tan(x)
    x_reduced = x;
    if (x_reduced < 0) {
        x_reduced = -x_reduced;
        symm_sign = -1;
    }
    // Reduce x to [0, 2pi)
    x_reduced = mpf_remainder(x_reduced, two_pi);
    // Furthur reduce x to [0, pi/2)
    if ((pi_over_2 < x_reduced) && (x_reduced <= _PI)) {
        x_reduced = _PI - x_reduced;
    }
    if ((_PI < x_reduced) && (x_reduced <= three * two_pi)) {
        x_reduced = three * two_pi - x_reduced;
        symm_sign *= -1;
    }
    if ((three * two_pi < x_reduced) && (x_reduced <= two_pi)) {
        x_reduced = two_pi - x_reduced;
        symm_sign *= -1;
    }
    // Calculate tan(x) using Taylor series
    mp_bitcnt_t k = std::floor(std::sqrt(req_precision / 2));
    // We need some additional precision for successive applications of the double-angle formula
    mp_bitcnt_t additional_precision = 0;
    if (k % 64 != 0) {
        additional_precision = ((k / 64) + 1) * 64;
    }
    mp_bitcnt_t _req_precision;
    _req_precision = req_precision + additional_precision;
    // Constants and variables
    mpf_class cosx(0.0, _req_precision);
    mpf_class sinx(0.0, _req_precision);
    mpf_class tanx(0.0, _req_precision);

    cosx = cos_taylor_reduced(x, true);
    sinx = sinx_from_cos_internal(x, true);
    tanx = (sinx / cosx);
    tanx *= symm_sign;
    return tanx;
}
inline mpf_class tan(const mpf_class &x) { return tan_from_sin_cos(x); }
inline mpf_class pow_from_exp_log(const mpf_class &x, const mpf_class &y) {
    mp_bitcnt_t req_precision = x.get_prec();
    mp_bitcnt_t req_precision_y = y.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
    assert(req_precision_y == mpf_get_default_prec());
#endif
    if (req_precision != req_precision_y) {
        throw std::runtime_error("Precision mismatch between x and y.");
    }
    mpf_class log_x(0.0, req_precision);
    mpf_class y_log_x(0.0, req_precision);
    mpf_class result(0.0, req_precision);
    log_x = log(mpf_class(x));
    y_log_x = y * log_x;
    result = exp(mpf_class(y_log_x));
    return result;
}
inline mpf_class pow(const mpf_class &x, const mpf_class &y) {
    if (mpf_integer_p(y.get_mpf_t()) != 0) {
        if (y > 0) {
            unsigned long int y_uint = mpf_get_ui(y.get_mpf_t());
            mpf_class result(0.0, x.get_prec());
            mpf_pow_ui(result.get_mpf_t(), x.get_mpf_t(), y_uint);
            return result;
        } else {
            // Handle y < 0 using exp and log
            unsigned long int y_uint = -mpf_get_ui(y.get_mpf_t());
            mpf_class result(0.0, x.get_prec());
            mpf_pow_ui(result.get_mpf_t(), x.get_mpf_t(), y_uint);
            return 1.0 / result;
        }
    } else {
        return pow_from_exp_log(x, y);
    }
}
inline mpf_class log2_from_log(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class _log2(0.0, req_precision);
    mpf_class result(0.0, req_precision);
    _log2 = log(mpf_class(2));
    result = log(x) / _log2;
    return result;
}
inline mpf_class log2(const mpf_class &x) { return log2_from_log(x); }
inline mpf_class log10_from_log(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class _log10(0.0, req_precision);
    mpf_class result(0.0, req_precision);
    _log10 = log(mpf_class(10));
    result = log(x) / _log10;
    return result;
}
inline mpf_class log10(const mpf_class &x) { return log10_from_log(x); }
inline mpf_class atan_AGM(const mpf_class &_x) {
    mp_bitcnt_t req_precision = _x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // Constants and variables
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class three(3.0, req_precision);
    mpf_class four(4.0, req_precision);
    mpf_class s(2.0, req_precision);
    mpf_class epsilon(2.0, req_precision);
    mpf_class x(zero, req_precision);
    mpf_class v(zero, req_precision);
    mpf_class q(one, req_precision);
    mpf_class ai(zero, req_precision);
    mpf_class bi(one, req_precision);
    mpf_class ci(zero, req_precision);
    mpf_class si(one, req_precision);
    mpf_class vi(one, req_precision);
    mpf_class qi(one, req_precision);
    mpf_class atanx(zero, req_precision);
    mpf_class atanx_refined(zero, req_precision);
    int sign = 1;
    int reduce = 1;
    if (_x < 0) {
        x = -_x;
        sign = -1;
    } else
        x = _x;
    //
    if (x >= 1) {
        x = (sqrt(one + x * x) - one) / x;
        reduce = 2;
    }
    s.div_2exp(req_precision / 2);
    epsilon.div_2exp(req_precision);

    v = x / (one + sqrt(one + x * x));
    q = one;

    si = s, vi = v, qi = q;
    while (one - si >= epsilon) {
        qi = two * qi / (one + si);
        ai = two * si * vi / (one + vi * vi);
        bi = ai / (one + sqrt(one - ai * ai));
        ci = (vi + bi) / (one - vi * bi);
        vi = ci / (one + sqrt(one + ci * ci));
        si = two * sqrt(si) / (one + si);
    }
    atanx = qi * log((one + vi) / (one - vi));
    atanx = atanx * sign * reduce;
    atanx_refined = atanx - (tan(atanx) - _x) / (one + tan(atanx) * tan(atanx));
    return atanx_refined;
}
inline mpf_class atan2(const mpf_class &y, const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
    mp_bitcnt_t req_precision_y = y.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    assert(req_precision == req_precision_y);
    mpf_class zero(0.0, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class three(3.0, req_precision);
    mpf_class four(4.0, req_precision);
    mpf_class pi_(0.0, req_precision);
    //    mpf_class pi_quarter_(0.0, req_precision);
    mpf_class pi_over_2_(0.0, req_precision);
    //    mpf_class pi_3_over_4_(0.0, req_precision);
    pi_ = const_pi(req_precision);
    //    pi_quarter_ = pi_ / four;
    pi_over_2_ = pi_ / two;
    //    pi_3_over_4_ = pi_ * three / four;
    if (y == 0 && x == 0) {
        return 0;
    }
    if (x == 0) {
        if (y > 0)
            return pi_over_2_;
        if (y < 0)
            return -pi_over_2_;
    }
    if (y == 0) {
        if (x > 0)
            return 0;
        if (x < 0)
            return pi_;
    }
    /* GMP does not comply IEEE 754 or IEC 60559, it doesn't have NaN nor +-Inf
    if (x == std::numeric_limits<mpf_class>::infinity()) {
        if (y == std::numeric_limits<mpf_class>::infinity())
            return pi_quarter_;
        if (y == -std::numeric_limits<mpf_class>::infinity())
            return -pi_quarter_;
        return zero;
    }
    if (x == -std::numeric_limits<mpf_class>::infinity()) {
        if (y == std::numeric_limits<mpf_class>::infinity())
            return pi_3_over_4_;
        if (y == -std::numeric_limits<mpf_class>::infinity())
            return -pi_3_over_4_;
        return pi_;
    }
    if (y == std::numeric_limits<mpf_class>::infinity())
        return pi_over_2_;
    if (y == -std::numeric_limits<mpf_class>::infinity())
        return -pi_over_2_;
    */
    if (x > 0) {
        return atan_AGM(y / x);
    } else if (x < zero && y >= zero) {
        return atan_AGM(y / x) + pi_;
    } else {
        return atan_AGM(y / x) - pi_;
    }
}
inline mpf_class atan(const mpf_class &x) { return atan_AGM(x); }
inline mpf_class asin_AGM(const mpf_class &x) {
    if (x < -1 || x > 1) {
        throw std::out_of_range("Error: x must be between -1 and 1.");
    }
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class one(1.0, req_precision);
    mpf_class half(0.5, req_precision);
    mpf_class pi_over_2(1.0, req_precision);
    mpf_class asin_refined(1.0, req_precision);
    mpf_class result(0.0, req_precision);
    pi_over_2 = const_pi(req_precision) * half;
    if (x == -1)
        return -pi_over_2;
    if (x == 1)
        return pi_over_2;
    mpf_class sqrt_one_minus_x2(0.0, req_precision);
    sqrt_one_minus_x2 = sqrt(one - x * x);
    result = atan_AGM(x / sqrt_one_minus_x2);
    asin_refined = result - (sin(result) - x) / cos(result);
    return asin_refined;
}

// arcsin(x) using Taylor series expansion around x = 0
// |x| \sim 1 there is a serious convergence problem
inline mpf_class arcsin_taylor(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class abs_x(0.0, req_precision);
    mpf_class term(0.0, req_precision);
    mpf_class arcsin_x(0.0, req_precision);
    mpf_class arcsin_x_prev(0.0, req_precision);
    mpf_class fact_2n(1.0, req_precision);
    mpf_class fact_n(1.0, req_precision);
    mpf_class pow_4n(1.0, req_precision);
    mpf_class two_n_plus_one(1.0, req_precision);
    mpf_class x_pow_2n_plus_1(0.0, req_precision);
    mpf_class epsilon(2.0, req_precision);
    mpf_class n(0.0, req_precision);

    mpf_class zero(0.0, req_precision);
    mpf_class quarter(0.25, req_precision);
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class four(4.0, req_precision);

    mpf_class pi_over_2(1.0, req_precision);
    mpf_class half(0.5, req_precision);
    pi_over_2 = const_pi(req_precision) * half;
    if (x == -1)
        return -pi_over_2;
    if (x == 1)
        return pi_over_2;

    bool negative = false;
    abs_x = x;
    if (x < 0) {
        negative = true;
        abs_x = -x;
    }
    term = abs_x;
    epsilon.div_2exp(req_precision);

    x_pow_2n_plus_1 = abs_x;                             // x^(2n+1)
    for (mp_bitcnt_t _n = 0; _n < req_precision; ++_n) { // |x| \sim 1 there is a serious convergence problem
        if (_n > 0) {
            n = mpf_class(_n);
            fact_2n *= (two * n) * (two * n - one);
            fact_n *= n;
            pow_4n *= four;
            two_n_plus_one = two * n + one;
            x_pow_2n_plus_1 *= abs_x * abs_x;
            term = (fact_2n / (pow_4n * fact_n * fact_n * two_n_plus_one)) * x_pow_2n_plus_1;
        }
        arcsin_x_prev = arcsin_x;
        arcsin_x += term;
        if (abs(arcsin_x - arcsin_x_prev) < epsilon) {
            break;
        }
    }
    return negative ? -arcsin_x : arcsin_x;
}
inline mpf_class asin(const mpf_class &x) { return asin_AGM(x); }
inline mpf_class acos(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class acosx(0.0, req_precision);
    mpf_class half_pi(0.0, req_precision);
    half_pi = const_pi(req_precision) * 0.5;
    acosx = half_pi - asin(x);
    return acosx;
}
inline mpf_class sinh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class exp_x(0.0, req_precision);
    mpf_class exp_neg_x(0.0, req_precision);
    exp_x = exp(x);
    exp_neg_x = exp(-x);
    return (exp_x - exp_neg_x) / 2;
}
inline mpf_class cosh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class exp_x(0.0, req_precision);
    mpf_class exp_neg_x(0.0, req_precision);
    exp_x = exp(x);
    exp_neg_x = exp(-x);
    return (exp_x + exp_neg_x) / 2;
}
inline mpf_class tanh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class sinh_x(0.0, req_precision);
    mpf_class cosh_x(0.0, req_precision);
    sinh_x = sinh(x);
    cosh_x = cosh(x);
    return sinh_x / cosh_x;
}
inline mpf_class asinh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    mpf_class one(1.0, req_precision);
    mpf_class result(0.0, req_precision);
    mpf_class sqrt_term(0.0, req_precision);
    sqrt_term = sqrt(x * x + one);
    result = log(x + sqrt_term);
    return result;
}
inline mpf_class acosh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    if (x < 1) {
        throw std::domain_error("acosh is not defined for x < 1");
    }
    mpf_class one(1.0, req_precision);
    mpf_class result(0.0, req_precision);
    mpf_class sqrt_term(0.0, req_precision);
    sqrt_term = sqrt(x * x - one);
    result = log(x + sqrt_term);
    return result;
}
inline mpf_class atanh(const mpf_class &x) {
    mp_bitcnt_t req_precision = x.get_prec();
#if defined ___GMPXX_MKII_NOPRECCHANGE___
    assert(req_precision == mpf_get_default_prec());
#endif
    // atanh is only defined for -1 < x < 1
    if (x <= -1 || x >= 1) {
        throw std::domain_error("atanh is not defined for |x| >= 1");
    }
    mpf_class one(1.0, req_precision);
    mpf_class two(2.0, req_precision);
    mpf_class result(0.0, req_precision);
    result = log((one + x) / (one - x)) / two;
    return result;
}
class gmp_randclass {
  public:
    // gmp_randinit_default, gmp_randinit_mt
    explicit gmp_randclass(void (*init_function)(gmp_randstate_t)) { init_function(state); }
    // gmp_randinit_lc_2exp
    gmp_randclass(void (*init_function)(gmp_randstate_t, mpz_srcptr, unsigned long, mp_bitcnt_t), mpz_class a, unsigned long c, mp_bitcnt_t m2exp) { init_function(state, a.get_mpz_t(), c, m2exp); }
    // gmp_randinit_lc_2exp_size
    gmp_randclass(int (*init_function)(gmp_randstate_t, mp_bitcnt_t), mp_bitcnt_t m2exp) {
        if (m2exp > 128) {
            throw std::length_error("m2exp parameter is too large");
        }
        if (init_function(state, m2exp) == 0) {
            throw std::runtime_error("Initialization failed");
        }
    }
    // gmp_randinit_set, only
    gmp_randclass(gmp_randalg_t alg, mp_bitcnt_t size) { __initialize(alg, size); }
    gmp_randclass(const gmp_randclass &) = delete;
    gmp_randclass &operator=(const gmp_randclass &) = delete;
    ~gmp_randclass() { gmp_randclear(state); }
    void seed(unsigned long int s) { gmp_randseed_ui(state, s); }
    void seed(const mpz_class &s) { gmp_randseed(state, s.get_mpz_t()); }
    mpz_class get_z_bits(mp_bitcnt_t bits) {
        mpz_class result;
        mpz_urandomb(result.get_mpz_t(), state, bits);
        return result;
    }
    mpz_class get_z_range(const mpz_class &n) {
        mpz_class result;
        mpz_urandomm(result.get_mpz_t(), state, n.get_mpz_t());
        return result;
    }
    mpf_class get_f(const mpf_class &op) {
        mpf_class result;
        mpf_urandomb(result.get_mpf_t(), state, op.get_prec());
        return result;
    }
    mpf_class get_f(mp_bitcnt_t prec = 0) {
#if defined ___GMPXX_MKII_NOPRECCHANGE___
        if (prec == 0)
            prec = mpf_get_default_prec();
#else
        if (prec == 0) {
            throw std::runtime_error("Prec must be specified.");
        }
#endif
        mpf_class result(0, prec);
        mpf_urandomb(result.get_mpf_t(), state, prec);
        return result;
    }

  private:
    gmp_randstate_t state;
    // obslated and equivalent to gmp_randinit_lc_2exp_size
    // cf. https://gmplib.org/manual/Random-State-Initialization
    void __initialize(gmp_randalg_t alg, mp_bitcnt_t size) {
        switch (alg) {
        case GMP_RAND_ALG_LC: // GMP_RAND_ALG_DEFAULT and 0 are the same as GMP_RAND_ALG_LC.
            if (gmp_randinit_lc_2exp_size(state, size) == 0) {
                throw std::runtime_error("Initialization failed");
            }
            break;
        default:
            gmp_errno |= GMP_ERROR_UNSUPPORTED_ARGUMENT;
            throw std::invalid_argument("Unsupported random algorithm");
        }
    }
}; // gmp_randclass
#if !defined ___GMPXX_DONT_USE_NAMESPACE___
} // namespace gmp
#endif

// mpf_class operator"" _mpf (const char *str)
// mpz_class operator"" _mpz (const char *str)
// mpq_class operator"" _mpq (const char *str)
#if defined ___GMPXX_DONT_USE_NAMESPACE___
inline mpz_class operator"" _mpz(unsigned long long int val) { return mpz_class(static_cast<unsigned long int>(val)); }
inline mpq_class operator"" _mpq(unsigned long long int val) { return mpq_class(static_cast<unsigned long int>(val), static_cast<unsigned long int>(1)); }
inline mpf_class operator"" _mpf(long double val) { return mpf_class(static_cast<double>(val)); }
#else
inline gmpxx::mpz_class operator"" _mpz(unsigned long long int val) { return gmpxx::mpz_class(static_cast<unsigned long int>(val)); }
inline gmpxx::mpq_class operator"" _mpq(unsigned long long int val) { return gmpxx::mpq_class(static_cast<unsigned long int>(val), static_cast<unsigned long int>(1)); }
inline gmpxx::mpf_class operator"" _mpf(long double val) { return gmpxx::mpf_class(static_cast<double>(val)); }
#endif
// in the manual, the following functions are avilable, but actually not.
// cf. https://gmplib.org/manual/C_002b_002b-Interface-Rationals
// "With C++11 compilers, integral rationals can be constructed with the syntax 123_mpq which is equivalent to mpq_class(123_mpz). Other rationals can be built as -1_mpq/2 or 0xb_mpq/123456_mpz."
#if !defined ___GMPXX_UDL_CHAR___
inline gmpxx::mpz_class operator"" _mpz(const char *str, [[maybe_unused]] std::size_t length) { return gmpxx::mpz_class(str); }
inline gmpxx::mpq_class operator"" _mpq(const char *str, [[maybe_unused]] std::size_t length) { return gmpxx::mpq_class(str); }
inline gmpxx::mpf_class operator"" _mpf(const char *str, [[maybe_unused]] std::size_t length) { return gmpxx::mpf_class(str); }
#endif

#if defined ___GMPXX_STRICT_COMPATIBILITY___
void swap(mpz_class &op1, mpz_class &op2) noexcept { op1.swap(op2); }
void swap(mpq_class &op1, mpq_class &op2) noexcept { op1.swap(op2); }
void swap(mpf_class &op1, mpf_class &op2) noexcept { op1.swap(op2); }
#endif

namespace std { // namespace std for numerical limits
template <>
#if defined ___GMPXX_DONT_USE_NAMESPACE___
class numeric_limits<mpz_class> {
#else
class numeric_limits<gmpxx::mpz_class> {
#endif
  public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr bool has_denorm_loss = false;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;
    static constexpr int digits = 0;
    static constexpr int digits10 = 0;
    static constexpr int max_digits10 = 0;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = round_toward_zero;
#if defined ___GMPXX_DONT_USE_NAMESPACE___
    static mpz_class min() noexcept { return mpz_class(0); }
    static mpz_class max() noexcept { return mpz_class(0); }
    static mpz_class lowest() noexcept { return mpz_class(0); }
    static mpz_class epsilon() noexcept { return mpz_class(0); }
    static mpz_class round_error() noexcept { return mpz_class(0); }
    static mpz_class infinity() noexcept { return mpz_class(0); }
    static mpz_class quiet_NaN() noexcept { return mpz_class(0); }
    static mpz_class signaling_NaN() noexcept { return mpz_class(0); }
    static mpz_class denorm_min() noexcept { return mpz_class(0); }
#else
    static gmpxx::mpz_class min() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class max() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class lowest() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class epsilon() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class round_error() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class infinity() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class quiet_NaN() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class signaling_NaN() noexcept { return gmpxx::mpz_class(0); }
    static gmpxx::mpz_class denorm_min() noexcept { return gmpxx::mpz_class(0); }
#endif
};
template <>
#if defined ___GMPXX_DONT_USE_NAMESPACE___
class numeric_limits<mpq_class> {
#else
class numeric_limits<gmpxx::mpq_class> {
#endif
  public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr bool has_denorm_loss = false;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;
    static constexpr int digits = 0;
    static constexpr int digits10 = 0;
    static constexpr int max_digits10 = 0;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = round_toward_zero;
#if defined ___GMPXX_DONT_USE_NAMESPACE___
    static mpq_class min() noexcept { return mpq_class(0); }
    static mpq_class max() noexcept { return mpq_class(0); }
    static mpq_class lowest() noexcept { return mpq_class(0); }
    static mpq_class epsilon() noexcept { return mpq_class(0); }
    static mpq_class round_error() noexcept { return mpq_class(0); }
    static mpq_class infinity() noexcept { return mpq_class(0); }
    static mpq_class quiet_NaN() noexcept { return mpq_class(0); }
    static mpq_class signaling_NaN() noexcept { return mpq_class(0); }
    static mpq_class denorm_min() noexcept { return mpq_class(0); }
#else
    static gmpxx::mpq_class min() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class max() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class lowest() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class epsilon() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class round_error() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class infinity() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class quiet_NaN() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class signaling_NaN() noexcept { return gmpxx::mpq_class(0); }
    static gmpxx::mpq_class denorm_min() noexcept { return gmpxx::mpq_class(0); }
#endif
};
template <>
#if defined ___GMPXX_DONT_USE_NAMESPACE___
class numeric_limits<mpf_class> {
#else
class numeric_limits<gmpxx::mpf_class> {
#endif
  public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;
    static constexpr int digits = 0;
    static constexpr int digits10 = 0;
    static constexpr int max_digits10 = 0;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = round_toward_zero;
#if defined ___GMPXX_DONT_USE_NAMESPACE___
    static mpf_class min() noexcept { return mpf_class(0); }
    static mpf_class max() noexcept { return mpf_class(0); }
    static mpf_class lowest() noexcept { return mpf_class(0); }
    static mpf_class epsilon() noexcept { return mpf_class(0); }
    static mpf_class round_error() noexcept { return mpf_class(0); }
    static mpf_class infinity() noexcept { return mpf_class(0); }
    static mpf_class quiet_NaN() noexcept { return mpf_class(0); }
    static mpf_class signaling_NaN() noexcept { return mpf_class(0); }
    static mpf_class denorm_min() noexcept { return mpf_class(0); }
#else
    static gmpxx::mpf_class min() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class max() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class lowest() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class epsilon() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class round_error() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class infinity() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class quiet_NaN() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class signaling_NaN() noexcept { return gmpxx::mpf_class(0); }
    static gmpxx::mpf_class denorm_min() noexcept { return gmpxx::mpf_class(0); }
#endif
};
} // namespace std

#endif // ___GMPXX_MKII_H___
