/*
 * polynomial.h
 *
 *  Created on: Jan 28, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "variable.h"
#include "variable_order.h"
#include "integer.h"
#include "polynomial_context.h"
#include "assignment.h"

/**
 * Polynomials incorporate the context and coefficient data. The also carry
 * flags so as not to re-do any computation.
 *
 * The external if the external flag is on during construction, the polynomial
 * is marked as external. This means that the context data will be attached and
 * before every operation, the polynomial will be reordered to in the right
 * order. The non-external polynomials are useful if you are doing lots of
 * intermediate computation to produce a final (external) polynomial.
 *
 * For non-external polynomials no checks are done during operations, so if the
 * order has changed the operations might fail or produce wrong results.
 *
 * In any case, the polynomial context should not be changed externally during
 * any operations.
 */
typedef struct {

  /** Construct a zero polynomial (does not attach the context) */
  void (*construct) (lp_polynomial_t* A, const lp_polynomial_context_t* ctx);

  /** Construct a simple polynomial c*x^n */
  void (*construct_simple) (lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const lp_integer_t* c, lp_variable_t x, unsigned n);

  /** Construct a copy of the given polynomial (does not attach the context). */
  void (*construct_copy) (lp_polynomial_t* A, const lp_polynomial_t* from);

  /** Destruct the polynomial. */
  void (*destruct) (lp_polynomial_t* A);

  /** Allocate a new polynomial (unconstructed) */
  lp_polynomial_t* (*alloc) (void);

  /** Allocate and construct a new polynomial */
  lp_polynomial_t* (*new) (const lp_polynomial_context_t* ctx);

  /** Makr the polynomial as external */
  void (*set_external) (lp_polynomial_t* A);

  /** Swap two polynomials. */
  void (*swap) (lp_polynomial_t* A1, lp_polynomial_t* A2);

  /** Assign the polynomial a given polynomial. */
  void (*assign) (lp_polynomial_t* A, const lp_polynomial_t* from);

  /** Returns the context of the polynomial */
  const lp_polynomial_context_t* (*context) (const lp_polynomial_t* A);

  /** Returns the degree of the polynomial (in it's top variable) */
  size_t (*degree) (const lp_polynomial_t* A);

  /** Returns the top variable of the polynomial */
  lp_variable_t (*top_variable) (const lp_polynomial_t* A);

  /** Puts the k-th coefficient of A into C */
  void (*get_coefficient) (lp_polynomial_t* C, const lp_polynomial_t* A, size_t k);

  /** Get the reductum of the polynomial (the polynomial withough the leading coefficient) */
  void (*reductum) (lp_polynomial_t* R, const lp_polynomial_t* A);

  /** Get the model-based reductum of the polynomial (the polynomial withough the leading coefficient) */
  void (*reductum_m) (lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m);

  /** Returns true if the polynomial is a constant */
  int (*is_constant) (const lp_polynomial_t* A);

  /** Returns true if the polynomial is 0 */
  int (*is_zero) (const lp_polynomial_t* A);

  /** returns the sign of the polynomial in the model */
  int (*sgn) (const lp_polynomial_t* A, const lp_assignment_t* m);

  /**
   * Compare the two polynomials in the ring. Not necessarily +/- 1, could be
   * any integer, only the sign matters.
   */
  int (*cmp) (const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Compare the type of two polynomials. Type if either a constant or a
   * polynomial, where comparison is done only by the lead variable comparison
   * and constants are smaller than other polynomials. This is not a total
   * ORDER.
   */
  int (*cmp_type) (const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Returns true if A1 divides A2. */
  int (*divides) (const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Prints the polynomial to the given stream. */
  int (*print) (const lp_polynomial_t* A, FILE* out);

  /** Returns the string representation of the polynomial. */
  char* (*to_string) (const lp_polynomial_t* A);

  /** Compute S = A1 + A2. */
  void (*add) (lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute S = A1 - A2 in the given ring. */
  void (*sub) (lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute N = -C. */
  void (*neg) (lp_polynomial_t* N, const lp_polynomial_t* A);

  /** Compute P = A1 * A2. */
  void (*mul) (lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Multiplication with x^n.
   */
  void (*shl) (lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n);

  /** Compute P = C^n. */
  void (*pow) (lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n);

  /** Compute S += A1*A2. */
  void (*add_mul) (lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute S -= A1*A2. */
  void (*sub_mul) (lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Reduce the polynomial A in Z[y,x] using B in Z[y,x] so that
   *
   *   P*A = Q*B + R
   *
   * and
   *
   *   P in Z[y]
   *   Q, R in Z[y,x]
   *
   * with
   *
   *   deg(R) < deg(B) or deg(R) = 0
   */
  void (*reduce) (const lp_polynomial_t* A, const lp_polynomial_t* B,
      lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R);

  /** Compute A1 = D*A2 (assumes that A2 divides A1). */
  void (*div) (lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute A1 = D*A2 + R (assumes that exact division). */
  void (*rem) (lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute a*A1 = D*A2 + R (pseudo remainder). */
  void (*prem) (lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute a*A1 = D*A2 + R (sparse pseudo remainder). */
  void (*sprem) (lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute A1 = D*A2 + R (assumes that exact division). */
  void (*divrem) (lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute A_d = A' (in the top variable). */
  void (*derivative) (lp_polynomial_t* A_d, const lp_polynomial_t* A);

  /** Compute the greatest common divisor gcd(A1, A2). */
  void (*gcd) (lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /** Compute the least common multiple lcm(A1, A2). */
  void (*lcm) (lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Compute the resultant of A1 and A2 in their top variable. Both A1 and A2
   * must be (non-trivial) polynomials over the same variable.
   */
  void (*resultant) (lp_polynomial_t* res, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Compute the principal subresultant coefficients (psc) of A1 and A1. Bot A1
   * and A2 must be (non-trivial) polynomials over the same variable. If
   *  deg(A1) = m, deg(A2) = n, and output will be of size min(m, n).
   */
  void (*psc) (lp_polynomial_t** psc, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

  /**
   * Get the square-free factorization of the given polynomial. It will allocate
   * the given arrays and return the polynomials and their multiplicities in
   * them.
   */
  void (*factor_square_free) (const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size);

} polynomial_ops_t;

extern const polynomial_ops_t polynomial_ops;
