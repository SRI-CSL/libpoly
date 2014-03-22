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

#include "upolynomial/upolynomial_internal.h"

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
  integer_t* coefficients;
} upolynomial_dense_t;

/**
 * Operations on dense polynomials.
 */
typedef struct {

  /** Construct a 0 polynomial of given capacity */
  void (*construct) (upolynomial_dense_t* p_d, size_t capacity);

  /**
   * Construct a dense polynomial from a given polynomial. The capacity is set
   * to the given capactity.
   */
  void (*construct_p) (upolynomial_dense_t* p_d, size_t capacity, const upolynomial_t* p);

  /**
   * Destructs the polynomial.
   */
  void (*destruct) (upolynomial_dense_t* p_d);

  /**
   * Swap two polynomials.
   */
  void (*swap) (upolynomial_dense_t* p_d, upolynomial_dense_t* q_d);

  /**
   * p_d = q_d
   */
  void (*assign) (upolynomial_dense_t* p_d, const upolynomial_dense_t* q_d);

  /**
   * p_d = 0
   */
  void (*clear) (upolynomial_dense_t* p_d);

  /**
   * Returns true if the polynomial is zero.
   */
  int (*is_zero) (const upolynomial_dense_t* p_d);

  /**
   * Returns the leading coefficient of p_d.
   */
  const integer_t* (*lead_coeff) (const upolynomial_dense_t* p_d);

  /**
   * Evaluate the polynomial at a rational point.
   */
  void (*evaluate_at_rational) (const upolynomial_dense_t* p_d, const rational_t* x, rational_t* value);

  /**
   * Evaluate the polynomial at a rational point.
   */
  void (*evaluate_at_dyadic_rational) (const upolynomial_dense_t* p_d, const dyadic_rational_t* x, dyadic_rational_t* value);

  /**
   * Returns the sign of the polynomial at the given rational point x. The
   * polynomial is assumed integer, otherwise the operation makes no sense.
   */
  int (*sgn_at_rational) (const upolynomial_dense_t* p_d, const rational_t* x);

  /**
   * Returns the sign of the polynomial at the given rational point x. The
   * polynomial is assumed integer, otherwise the operation makes no sense.
   */
  int (*sgn_at_dyadic_rational) (const upolynomial_dense_t* p_d, const dyadic_rational_t* x);

  /**
   * Returns the sign of the polynomial at +inf.
   */
  int (*sgn_at_plus_inf) (const upolynomial_dense_t* p_d);

  /**
   * Returns the sign of the polynomial at -inf.
   */
  int (*sgn_at_minus_inf) (const upolynomial_dense_t* p_d);

  /**
   * Returns the sparse polynomial.
   */
  upolynomial_t* (*to_upolynomial)(const upolynomial_dense_t* p_d, int_ring K);

  /**
   * Print the polynomial.
   */
  int (*print) (const upolynomial_dense_t* p_d, FILE* file);

  /**
   * Call when modifying a coefficient, so as to keep internal consistency.
   */
  void (*touch) (upolynomial_dense_t* p_d, int_ring K, size_t degree);

  /**
   * Make the polynomial primitive (divide by gcd). If the positive flag is
   * true, the leading coefficient is positive will be made positive.
   */
  void (*mk_primitive_Z) (upolynomial_dense_t* p_d, int positive);

  /**
   * p_d *= c
   */
  void (*mult_c) (upolynomial_dense_t* p_d, int_ring K, const integer_t* c);

  /**
   * p_d /= c (c should should divide all coefficiants).
   */
  void (*div_c) (upolynomial_dense_t* p_d, int_ring K, const integer_t* c);

  /**
   * p_d += p*c
   */
  void (*add_mult_p_c) (upolynomial_dense_t* p_d, const upolynomial_t* p, const integer_t* c);

  /**
   * p_d += p*c
   */
  void (*add_mult_p_int) (upolynomial_dense_t* p_d, const upolynomial_t* p, int c);

  /**
   * p_d += p*m
   */
  void (*add_mult_p_mon) (upolynomial_dense_t* p_d, const upolynomial_t* p, const umonomial_t* m);

  /**
   * p_d -= p*m
   */
  void (*sub_mult_p_mon) (upolynomial_dense_t* p_d, const upolynomial_t* p, const umonomial_t* m);

  /**
   * p_d -= p*m
   */
  void (*sub_mult_mon) (upolynomial_dense_t* p_d, int_ring K, const upolynomial_dense_t* p, const umonomial_t* m);

  /**
   * p_d = -p_d
   */
  void (*negate) (upolynomial_dense_t* p_d, int_ring K);

  /**
   * p_d -= p*q
   */
  void (*sub_mult) (upolynomial_dense_t* p_d, int_ring K, const upolynomial_dense_t* p, const upolynomial_dense_t* q);

  /**
   * General division p = div*q + rem in the ring K, if exact. If not exact, then
   * we compute lcm(q)^(p_deg - q_deg + 1) p = div*q + rem.
   */
  void (*div_general) (int_ring K, int exact, const upolynomial_dense_t* p, const upolynomial_dense_t* q, upolynomial_dense_t* div, upolynomial_dense_t* rem);

  /**
   * Reduce a polynomial p using q in Z[x]. Result is
   *    a*p = div*q + red
   */
  void (*reduce_Z) (const upolynomial_dense_t* p_d, const upolynomial_dense_t* q_d, integer_t* a, upolynomial_dense_t* red_d);

  /**
   * Derivative.
   */
  void (*derivative) (int_ring K, const upolynomial_dense_t* p_d, upolynomial_dense_t* p_d_prime);

} upolynomial_dense_ops_t;

/** Implementation of the operations */
extern const upolynomial_dense_ops_t upolynomial_dense_ops;
