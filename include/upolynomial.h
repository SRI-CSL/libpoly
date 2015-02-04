/*
 * univariate_polynomial.h
 *
 *  Created on: Oct 28, 2013
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "integer.h"
#include "rational.h"
#include "interval.h"

/** Interface of the polynomial operations */
typedef struct {

  /**
   * Construct the polynomial given its coefficients. Coefficients should be
   * indexed by degree and they will be normalized according to the given ring.
   */
  lp_upolynomial_t* (*construct) (lp_int_ring K, size_t degree, const lp_integer_t* coefficients);

  /**
   * Construct the polynomial c*x^d.
   */
  lp_upolynomial_t* (*construct_power) (lp_int_ring K, size_t degree, long c);

  /**
   * Construct the polynomial given its coefficients. Coefficients should be
   * indexed by degree and they will be normalize according to the given ring.
   */
  lp_upolynomial_t* (*construct_from_int) (lp_int_ring K, size_t degree, const int* coefficients);

  /**
   * Construct the polynomial given its coefficients. Coefficients should be
   * indexed by degree and they will be normalize according to the given ring.
   */
  lp_upolynomial_t* (*construct_from_long) (lp_int_ring K, size_t degree, const long* coefficients);

  /**
   * Construct a copy of the polynomial.
   */
  lp_upolynomial_t* (*construct_copy) (const lp_upolynomial_t* p);

  /**
   * Construct a copy of the polynomial, but change the ring.
   */
  lp_upolynomial_t* (*construct_copy_K) (lp_int_ring K, const lp_upolynomial_t* p);

  /**
   * Frees the polynomial data and detaches the ring. */
  void (*delete) (lp_upolynomial_t* p);

  /**
   * Returns the degree of the polynomial. Note that the degree of the constat
   * 0 is 0.
   */
  size_t (*degree) (const lp_upolynomial_t* p);

  /**
   * Returns the field of the polynomial (unatached).
   */
  lp_int_ring (*ring) (const lp_upolynomial_t* p);

  /**
   * Sets the ring to given ring (has to be "larger" than existing).
   */
  void (*set_ring) (lp_upolynomial_t* p, lp_int_ring K);

  /**
   * Returns the lead coefficient of the given polynomial.
   */
  const lp_integer_t* (*lead_coeff) (const lp_upolynomial_t* p);

  /**
   * Unpack the polynomial into a dense representation. The out vector is
   * assumed to be large enough, and filled with 0 (only non-zero coefficients
   * will be copied into out).
   */
  void (*unpack) (const lp_upolynomial_t* p, lp_integer_t* out);

  /**
   * Print the polynomial to the output stream.
   */
  int (*print) (const lp_upolynomial_t* p, FILE* out);

  /**
   * Get a string representation of the polynomial (you own the memory).
   */
  char* (*to_string) (const lp_upolynomial_t* p);

  /**
   * Returns true if this is a zero polynomial
   */
  int (*is_zero) (const lp_upolynomial_t* p);

  /**
   * Returns true if this is the polynomial 1.
   */
  int (*is_one) (const lp_upolynomial_t* p);

  /**
   * Returns true if the polynomial is monic.
   */
  int (*is_monic) (const lp_upolynomial_t* p);

  /**
   * Returns true if the polynomial is primitive, i.e. gcd of all coefficients
   * is 1 and the leading coefficient is positive.
   */
  int (*is_primitive) (const lp_upolynomial_t* p);

  /**
   * Evaluates the polynomial.
   */
  void (*evaluate_at_integer) (const lp_upolynomial_t* p, const lp_integer_t* x, lp_integer_t* value);

  /**
   * Evaluates the polynomial. Only makes sense for polynomials in Z[x].
   */
  void (*evaluate_at_rational) (const lp_upolynomial_t* p, const lp_rational_t* x, lp_rational_t* value);

  /**
   * Evaluates the polynomial. Only makes sense for polynomials in Z[x].
   */
  void (*evaluate_at_dyadic_rational) (const lp_upolynomial_t* p, const lp_dyadic_rational_t* x, lp_dyadic_rational_t* value);

  /**
   * Get the sign of the polynomial in the given integer point.
   */
  int (*sgn_at_integer) (const lp_upolynomial_t* p, const lp_integer_t* x);

  /**
   * Get the sign of the polynomial in the given rational point. Only makes sense
   * for polynomials in Z[x].
   */
  int (*sgn_at_rational) (const lp_upolynomial_t* p, const lp_rational_t* x);

  /**
   * Get the sign of the polynomial in the given rational point. Only makes sense
   * for polynomials in Z[x].
   */
  int (*sgn_at_dyadic_rational) (const lp_upolynomial_t* p, const lp_dyadic_rational_t* x);

  /**
   * Compares two polynomials (lexicographic from highest coefficient) and
   * returns -1 if p < q, 0 if p == q, and 1 if p > q.
   */
  int (*cmp) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Add two polynomials (all operations in the same ring).
   */
  lp_upolynomial_t* (*add) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Subtract two polynomials (all operations in the same ring).
   */
  lp_upolynomial_t* (*sub) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Multiply the two polynomials (all operations in the same ring).
   */
  lp_upolynomial_t* (*multiply) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Multiply the two polynomials (all operations in the ring of p).
   */
  lp_upolynomial_t* (*multiply_c) (const lp_upolynomial_t* p, const lp_integer_t* c);

  /**
   * Power of a polynomial (all operations in the ring of p).
   */
  lp_upolynomial_t* (*power) (const lp_upolynomial_t* p, long pow);

  /**
   * Returns the derivative of the given polynomial. Note that deg(p') can be
   * less than deg(p) - 1 in some rings. For example in Z_4 (2x^2)' = 0.
   */
  lp_upolynomial_t* (*derivative) (const lp_upolynomial_t* p);

  /***
   * Returns true if p divides q.
   */
  int (*divides) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Return a polynomial with all the monomial degrees divided by the given
   * positive number a. All degrees must be divisible by a.
   */
  lp_upolynomial_t* (*div_degrees) (const lp_upolynomial_t* p, size_t a);

  /**
   * Returns the exact division of two polynomials. This assumes that p and q
   * are in the same ring and all needed coefficient inverses can be computed.
   */
  lp_upolynomial_t* (*div_exact) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Returns the exact division of the polynomial with a constant c. This
   * assumes that all coefficients of p are divisible by c.
   */
  lp_upolynomial_t* (*div_exact_c) (const lp_upolynomial_t* p, const lp_integer_t* c);

  /**
   * Returns the exact remainder of two polynomials. This assumes that p and q
   * are in the same ring all needed coefficient inverses can be computed.
   */
  lp_upolynomial_t* (*rem_exact) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Returns the exact division of two polynomials. This assumes that p and q
   * are in the same ring and all needed coefficient inverses can be computed.
   * The output in div and rem will be newly allocated.
   */
  void (*div_rem_exact) (const lp_upolynomial_t* p, const lp_upolynomial_t* q,
      lp_upolynomial_t** div, lp_upolynomial_t** rem);

  /**
   * Psuedo-division of polynomials, div and rem such that
   *
   *   lcm(q)^(p_deg - q_deg + 1) p = div*q + rem
   *
   * This assynes that deg(p) >= deg(q).
   *
   * Note: all computation is done in ring of p and q, but the algorithm doesn't
   * take advantage of existence of possible inverses -- algorithm proceeds as
   * if done in Z with individual operations performed in the ring.
   */
  void (*div_pseudo) (lp_upolynomial_t** div, lp_upolynomial_t** rem, const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Compute the content of the polynomial. Content of the polynomial is the
   * gcd of the coefficients, of the same sign as the leading coefficient.
   * NOTE: The gcd is computed in Z so p must be in p
   */
  void (*content_Z) (const lp_upolynomial_t* p, lp_integer_t* content);

  /**
   * Make the polynomial primitive. The polynomial is primitive if the content
   * is 1.
   */
  void (*make_primitive_Z) (lp_upolynomial_t* p);

  /**
   * Get the primitive part of the polynomial, the primitive part is p/content
   * and always has leading coefficient > 0.
   */
  lp_upolynomial_t* (*primitive_part) (const lp_upolynomial_t* p);

  /**
   * Computes the polynomial greatest common divisor of p and q. The rings of p
   * and q are assumed to be the same. The lc(gcd) > 0, and if the ring of p
   * and q is a field, it will be normalized to monic -- lc(gcd) == 1.
   */
  lp_upolynomial_t* (*gcd) (const lp_upolynomial_t* p, const lp_upolynomial_t* q);

  /**
   * Computes the extended gcd of p and q. The rings of p
   * and q are assumed to be the same. The lc(gcd) > 0, and if the ring of p
   * and q is a field, it will be normalized to monic -- lc(gcd) == 1.
   */
  lp_upolynomial_t* (*extended_gcd) (const lp_upolynomial_t* p, const lp_upolynomial_t* q, lp_upolynomial_t** u, lp_upolynomial_t** v);

  /**
   * Given p, q, and r solve the equation
   *
   *  u*p+v*q = r
   *
   * for u an d. Assumes that gcd(p, q) divides r. Result such that
   *
   *   deg(u) < deg(q), deg(v) < deg(p)
   */
  void (*solve_bezout) (const lp_upolynomial_t* p, const lp_upolynomial_t* q, const lp_upolynomial_t* r,
      lp_upolynomial_t** u, lp_upolynomial_t** v);

  /**
   * Returns the factorization of the given polynomial in its ring.
   */
  lp_upolynomial_factors_t* (*factor) (const lp_upolynomial_t* p);

  /**
   * Returns the square-free factorization of the given polynomial in its ring.
   * In a square-free factorization each factor is square-free. Individual
   * factors are also mutually prime, i.e. gcd(f_i, f_j) = 1 for i != j.
   */
  lp_upolynomial_factors_t* (*factor_square_free) (const lp_upolynomial_t* p);

  /**
   * Return the Sturm sequence of the given polynomial. The arrays S will be
   * allocated, and the user should de-allocate it. The size parameter will be
   * updated with the size of the array.
   */
  void (*sturm_sequence) (const lp_upolynomial_t* f, lp_upolynomial_t*** S, size_t* size);

  /**
   * Counts the number of real roots in the given interval. If the interval is
   * 0, it counts through (-inf, inf).
   */
  int (*roots_count) (const lp_upolynomial_t* p, const interval_t* ab);

  /**
   * Isolate the distinct real roots of the given polynomial.
   */
  void (*roots_isolate) (const lp_upolynomial_t* p, lp_algebraic_number_t* roots, size_t* roots_size);

} upolynomial_ops_struct;

/** Function table of the polynomial operations */
extern const upolynomial_ops_struct upolynomial_ops;

/** Interface to the polynomial factors structure */
typedef struct {

  /** Construct the factors with no factors, and constant 1. */
  lp_upolynomial_factors_t* (*construct) (void);

  /**
   * Free the memory and the factor polynomials. If destruct_factors is true
   * then the individual factors are also destructed (should be, unless you
   * copied the factors somewhere else).
   */
  void (*destruct) (lp_upolynomial_factors_t* f, int destruct_factors);

  /** Clear the factor polynomials */
  void (*clear) (lp_upolynomial_factors_t* f);

  /** Swap the two factorizations */
  void (*swap) (lp_upolynomial_factors_t* f1, lp_upolynomial_factors_t* f2);

  /** Get the number of factors */
  size_t (*size) (const lp_upolynomial_factors_t* f);

  /** Get a factor with the given index i < size() */
  lp_upolynomial_t* (*get_factor) (lp_upolynomial_factors_t* f, size_t i, size_t* multiplicity);

  /** Returns the constant of the factorization */
  const lp_integer_t* (*get_constant) (const lp_upolynomial_factors_t* f);

  /** Add a factor with the given degree */
  void (*add) (lp_upolynomial_factors_t* f, lp_upolynomial_t* p, size_t d);

  /** Print the factors */
  int (*print) (const lp_upolynomial_factors_t* f, FILE* out);

  /** Get the ring */
  lp_int_ring (*ring) (const lp_upolynomial_factors_t* f);

  /**
   * Set the ring of all polynomials to K. This is only possible if K is
   * "larger" than the existing ring.
   */
  void (*set_ring) (lp_upolynomial_factors_t* f, lp_int_ring K);

} upolynomial_factors_ops_t;

extern const upolynomial_factors_ops_t upolynomial_factors_ops;

