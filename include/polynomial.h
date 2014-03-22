/*
 * polynomial.h
 *
 *  Created on: Jan 28, 2014
 *      Author: dejan
 */

#pragma once

#include "variable.h"
#include "variable_order.h"
#include "integer.h"
#include "polynomial_context.h"
#include "assignment.h"

typedef struct polynomial_rec_struct polynomial_rec_t;

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
typedef struct polynomial_struct polynomial_t;

typedef struct {

  /** Construct a zero polynomial (does not attach the context) */
  void (*construct) (polynomial_t* A, const polynomial_context_t* ctx, int external);

  /** Construct a simple polynomial c*x^n */
  void (*construct_simple) (polynomial_t* A, const polynomial_context_t* ctx, int external, const integer_t* c, variable_t x, unsigned n);

  /** Construct a copy of the given polynomial (does not attach the context). */
  void (*construct_copy) (polynomial_t* A, const polynomial_t* from, int external);

  /** Destruct the polynomial. */
  void (*destruct) (polynomial_t* A);

  /** Allocate a new polynomial (unconstructed) */
  polynomial_t* (*alloc) (void);

  /** Allocate and construct a new polynomial */
  polynomial_t* (*new) (const polynomial_context_t* ctx, int external);

  /** Swap two polynomials. */
  void (*swap) (polynomial_t* A1, polynomial_t* A2);

  /** Assign the polynomial a given polynomial. */
  void (*assign) (polynomial_t* A, const polynomial_t* from);

  /** Returns the context of the polynomial */
  const polynomial_context_t* (*context) (const polynomial_t* A);

  /** Returns the degree of the polynomial (in it's top variable) */
  size_t (*degree) (const polynomial_t* A);

  /** Returns the top variable of the polynomial */
  variable_t (*top_variable) (const polynomial_t* A);

  /** Puts the k-th coefficient of A into C */
  void (*get_coefficient) (polynomial_t* C, const polynomial_t* A, size_t k);

  /** Returns true if the polynomial is a constant */
  int (*is_constant) (const polynomial_t* A);

  /** Returns true if the polynomial is 0 */
  int (*is_zero) (const polynomial_t* A);

  /** returns the sign of the polynomial in the model */
  int (*sgn) (const polynomial_t* A, const assignment_t* m);

  /**
   * Compare the two polynomials in the ring. Not necessarily +/- 1, could be
   * any integer, only the sign matters.
   */
  int (*cmp) (const polynomial_t* A1, const polynomial_t* A2);

  /**
   * Compare the type of two polynomials. Type if either a constant or a
   * polynomial, where comparison is done only by the lead variable comparison
   * and constants are smaller than other polynomials. This is not a total
   * ORDER.
   */
  int (*cmp_type) (const polynomial_t* A1, const polynomial_t* A2);

  /** Returns true if A1 divides A2. */
  int (*divides) (const polynomial_t* A1, const polynomial_t* A2);

  /** Prints the polynomial to the given stream. */
  int (*print) (const polynomial_t* A, FILE* out);

  /** Returns the string representation of the polynomial. */
  char* (*to_string) (const polynomial_t* A);

  /** Compute S = A1 + A2. */
  void (*add) (polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute S = A1 - A2 in the given ring. */
  void (*sub) (polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute N = -C. */
  void (*neg) (polynomial_t* N, const polynomial_t* A);

  /** Compute P = A1 * A2. */
  void (*mul) (polynomial_t* P, const polynomial_t* A1, const polynomial_t* A2);

  /**
   * Multiplication with x^n.
   */
  void (*shl) (polynomial_t* S, const polynomial_t* A, unsigned n);

  /** Compute P = C^n. */
  void (*pow) (polynomial_t* P, const polynomial_t* A, unsigned n);

  /** Compute S += A1*A2. */
  void (*add_mul) (polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute S -= A1*A2. */
  void (*sub_mul) (polynomial_t* S, const polynomial_t* A1, const polynomial_t* A2);

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
  void (*reduce) (const polynomial_t* A, const polynomial_t* B,
      polynomial_t* P, polynomial_t* Q, polynomial_t* R);

  /** Compute A1 = D*A2 (assumes that A2 divides A1). */
  void (*div) (polynomial_t* D, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute A1 = D*A2 + R (assumes that exact division). */
  void (*rem) (polynomial_t* R, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute A1 = D*A2 + R (assumes that exact division). */
  void (*divrem) (polynomial_t* D, polynomial_t* R, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute A_d = A' (in the top variable). */
  void (*derivative) (polynomial_t* A_d, const polynomial_t* A);

  /** Compute the greatest common divisor gcd(A1, A2). */
  void (*gcd) (polynomial_t* gcd, const polynomial_t* A1, const polynomial_t* A2);

  /** Compute the least common multiple lcm(A1, A2). */
  void (*lcm) (polynomial_t* lcm, const polynomial_t* A1, const polynomial_t* A2);

  /**
   * Compute the resultant of A1 and A2 in their top variable. Both A1 and A2
   * must be (non-trivial) polynomials over the same variable.
   */
  void (*resultant) (polynomial_t* res, const polynomial_t* A1, const polynomial_t* A2);

  /**
   * Compute the principal subresultant coefficients (psc) of A1 and A1. Bot A1
   * and A2 must be (non-trivial) polynomials over the same variable. If
   *  deg(A1) = m, deg(A2) = n, and output will be of size min(m, n).
   */
  void (*psc) (polynomial_t** psc, const polynomial_t* A1, const polynomial_t* A2);

  /** Set the power symbol for printouts */
  void (*set_power_symbol) (const char* power);

} polynomial_ops_t;

extern const polynomial_ops_t polynomial_ops;
