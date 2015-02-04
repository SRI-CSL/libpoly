/*
 * upolynomial_buffer.h
 *
 *  Created on: Nov 17, 2013
 *      Author: dejan
 */

#pragma once

#include <stddef.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"
#include "upolynomial/umonomial.h"

#include <upolynomial.h>

/**
 * Buffer for dense polynomial operations. When computing with a buffer the
 * ring to compute in is always passed in (unless one of the arguments has a
 * ring associated with it). If you change the coefficients manually notify
 * it with a call to touch. The capacity of the polynomial is never changed so
 * make sure you allocate enough.
 */
typedef struct {
  /** Current capacity of the buffer */
  size_t capacity;
  /** Current size of the polynomial (first non-zero coefficient, unless 0, then 1) */
  size_t size;
  /** The array of coefficients */
  lp_integer_t* coefficients;
} upolynomial_dense_t;

/** Construct a 0 polynomial of given capacity */
void upolynomial_dense_construct(upolynomial_dense_t* p_d, size_t capacity);

/**
 * Construct a dense polynomial from a given polynomial. The capacity is set
 * to the given capactity.
 */
void upolynomial_dense_construct_p(upolynomial_dense_t* p_d, size_t capacity,
    const lp_upolynomial_t* p);

/**
 * Destructs the polynomial.
 */
void upolynomial_dense_destruct(upolynomial_dense_t* p_d);

/**
 * Swap two polynomials.
 */
void upolynomial_dense_swap(upolynomial_dense_t* p_d, upolynomial_dense_t* q_d);

/**
 * p_d = q_d
 */
void upolynomial_dense_assign(upolynomial_dense_t* p_d, const upolynomial_dense_t* q_d);

/**
 * p_d = 0
 */
void upolynomial_dense_clear(upolynomial_dense_t* p_d);

/**
 * Returns true if the polynomial is zero.
 */
int upolynomial_dense_is_zero(const upolynomial_dense_t* p_d);

/**
 * Returns the leading coefficient of p_d.
 */
const lp_integer_t* upolynomial_dense_lead_coeff(const upolynomial_dense_t* p_d);

/**
 * Evaluate the polynomial at a rational point.
 */
void upolynomial_dense_evaluate_at_rational(const upolynomial_dense_t* p_d, const lp_rational_t* x,
    lp_rational_t* value);

/**
 * Evaluate the polynomial at a rational point.
 */
void upolynomial_dense_evaluate_at_dyadic_rational(const upolynomial_dense_t* p_d,
    const lp_dyadic_rational_t* x, lp_dyadic_rational_t* value);

/**
 * Returns the sign of the polynomial at the given rational point x. The
 * polynomial is assumed integer, otherwise the operation makes no sense.
 */
int upolynomial_dense_sgn_at_rational(const upolynomial_dense_t* p_d, const lp_rational_t* x);

/**
 * Returns the sign of the polynomial at the given rational point x. The
 * polynomial is assumed integer, otherwise the operation makes no sense.
 */
int upolynomial_dense_sgn_at_dyadic_rational(const upolynomial_dense_t* p_d,
    const lp_dyadic_rational_t* x);

/**
 * Returns the sign of the polynomial at +inf.
 */
int upolynomial_dense_sgn_at_plus_inf(const upolynomial_dense_t* p_d);

/**
 * Returns the sign of the polynomial at -inf.
 */
int upolynomial_dense_sgn_at_minus_inf(const upolynomial_dense_t* p_d);

/**
 * Returns the sparse polynomial.
 */
lp_upolynomial_t* upolynomial_dense_to_upolynomial(const upolynomial_dense_t* p_d, lp_int_ring K);

/**
 * Call when modifying a coefficient, so as to keep internal consistency.
 */
void upolynomial_dense_touch(upolynomial_dense_t* p_d, size_t degree);

/**
 * Make the polynomial primitive (divide by gcd). If the positive flag is
 * true, the leading coefficient is positive will be made positive.
 */
void upolynomial_dense_mk_primitive_Z(upolynomial_dense_t* p_d, int positive);

/**
 * p_d *= c
 */
void upolynomial_dense_mult_c(upolynomial_dense_t* p_d, lp_int_ring K, const lp_integer_t* c);

/**
 * p_d /= c (c should should divide all coefficiants).
 */
void upolynomial_dense_div_c(upolynomial_dense_t* p_d, lp_int_ring K, const lp_integer_t* c);

/**
 * p_d += p*c
 */
void upolynomial_dense_add_mult_p_c(upolynomial_dense_t* p_d, const lp_upolynomial_t* p,
    const lp_integer_t* c);

/**
 * p_d += p*c
 */
void upolynomial_dense_add_mult_p_int(upolynomial_dense_t* p_d, const lp_upolynomial_t* p, int c);

/**
 * p_d += p*m
 */
void upolynomial_dense_add_mult_p_mon(upolynomial_dense_t* p_d, const lp_upolynomial_t* p,
    const umonomial_t* m);

/**
 * p_d -= p*m
 */
void upolynomial_dense_sub_mult_p_mon(upolynomial_dense_t* p_d, const lp_upolynomial_t* p,
    const umonomial_t* m);

/**
 * p_d -= p*m
 */
void upolynomial_dense_sub_mult_mon(upolynomial_dense_t* p_d, lp_int_ring K,
    const upolynomial_dense_t* p, const umonomial_t* m);

/**
 * p_d = -p_d
 */
void upolynomial_dense_negate(upolynomial_dense_t* p_d, lp_int_ring K);

/**
 * p_d -= p*q
 */
void upolynomial_dense_sub_mult(upolynomial_dense_t* p_d, lp_int_ring K,
    const upolynomial_dense_t* p, const upolynomial_dense_t* q);

/**
 * General division p = div*q + rem in the ring K, if exact. If not exact, then
 * we compute lcm(q)^(p_deg - q_deg + 1) p = div*q + rem.
 */
void upolynomial_dense_div_general(lp_int_ring K, int exact, const upolynomial_dense_t* p,
    const upolynomial_dense_t* q, upolynomial_dense_t* div,
    upolynomial_dense_t* rem);

/**
 * Reduce a polynomial p using q in Z[x]. Result is
 *    a*p = div*q + red
 */
void upolynomial_dense_reduce_Z(const upolynomial_dense_t* p_d, const upolynomial_dense_t* q_d,
    lp_integer_t* a, upolynomial_dense_t* red_d);

/**
 * Derivative.
 */
void upolynomial_dense_derivative(lp_int_ring K, const upolynomial_dense_t* p_d,
    upolynomial_dense_t* p_d_prime);

