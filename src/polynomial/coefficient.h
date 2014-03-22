/*
 * polynomial_coefficient.h
 *
 *  Created on: Feb 4, 2014
 *      Author: dejan
 */

#pragma once

#include "variable.h"
#include "assignment.h"

#include "number/integer.h"

/** Type of the coefficients */
typedef enum {
  /** A purely numeric coefficient */
  COEFFICIENT_NUMERIC,
  /** A recursive polynomial coefficient */
  COEFFICIENT_POLYNOMIAL
} coefficient_type_t;

typedef struct polynomial_rec_struct polnomial_rec_t;
typedef struct coefficient_struct coefficient_t;

/** Recursive nodes in the tree representation of the polynomial */
typedef struct polynomial_rec_struct {
  /** The used size of the coefficient array */
  size_t size     : 16;
  /** Capacity of the coefficient array */
  size_t capacity : 16;
  /** The main variable */
  variable_t x    : 32;
  /** Coefficients */
  coefficient_t* coefficients;
} polynomial_rec_t;

/**
 * Value of each coefficient is either a base value or a polynomial.
 */
typedef union {
  integer_t num;
  polynomial_rec_t rec;
} coefficient_union_t;

/**
 * A coefficient is a tagged union of coefficient values.
 */
typedef struct coefficient_struct {
  coefficient_type_t type;
  coefficient_union_t value;
} coefficient_t;

#define SIZE(C) ((C)->value.rec.size)
#define CAPACITY(C) ((C)->value.rec.capacity)
#define COEFF(C, i) ((C)->value.rec.coefficients + (i))
#define VAR(C) ((C)->value.rec.x)

/** Context for the polynomial operations */
typedef struct polynomial_context_struct polynomial_context_t;

/**
 * Type of remaindering in the reduce method.
 */
typedef enum {
    /** The standard dense pseudo-remainder */
    REMAINDERING_PSEUDO_DENSE,
    /** Assume the leading coefficients are divisible by the lc of the divisor */
    REMAINDERING_EXACT_SPARSE,
    /** The sparse pseudo-remainder */
    REMAINDERING_PSEUDO_SPARSE,
    /** Use the LCM of the leading coefficients */
    REMAINDERING_LCM_SPARSE
} remaindering_type_t;

/**
 * Interface to operations on coefficients. All operations can take a polynomial
 * context and can take the same coefficient as input output (i.e. it is
 * permissible to write add(ctx, a, a, a)).
 *
 * Most operations, if not indicated differently, assume that the the input
 * coefficients (const) all share the same context, i.e. same base ring and the same
 * order of variables.
 */
typedef struct {

  /** Construct a zero coefficient */
  void (*construct) (const polynomial_context_t* ctx, coefficient_t* C);

  /** Construct a coefficient from integer */
  void (*construct_from_int) (const polynomial_context_t* ctx, coefficient_t* C, long C_int);

  /** Construct a coefficient from integer */
  void (*construct_from_integer) (const polynomial_context_t* ctx, coefficient_t* C, const integer_t* C_integer);

  /** Construct from a univariate polynomial */
  void (*construct_from_univariate) (const polynomial_context_t* ctx, coefficient_t* C, const upolynomial_t* p, variable_t x);

  /** Construct a simple monomial coefficient a*x^n */
  void (*construct_simple) (const polynomial_context_t* ctx, coefficient_t* C, const integer_t* a, variable_t x, unsigned n);

  /** Construct a copy of the given coefficient. */
  void (*construct_copy) (const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

  /** Destructs the coefficient. */
  void (*destruct) (coefficient_t* C);

  /** Swap two coefficients. */
  void (*swap) (coefficient_t* C1, coefficient_t* C2);

  /** Assign the coefficient a given coefficient. */
  void (*assign) (const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

  /** Assign the coefficient a give integer */
  void (*assign_int) (const polynomial_context_t* ctx, coefficient_t* C, long x);

  /** Returns the univariate version of this coefficient (must be univariate) */
  upolynomial_t* (*to_univariate) (const polynomial_context_t* ctx, const coefficient_t* C);

  /** Returns true if the coefficient is a constant */
  int (*is_constant) (const coefficient_t* C);

  /** The degree of the coefficient */
  size_t (*degree) (const coefficient_t* C);

  /** The top variable of the coefficient */
  variable_t (*top_variable) (const coefficient_t* C);

  /** Get the k-th coefficient of this coefficient */
  const coefficient_t* (*get_coefficient) (const coefficient_t* C, size_t k);

  /** Get the leading coefficient of this coefficient */
  const coefficient_t* (*lc) (const coefficient_t* C);

  /** Returns true if the coefficient is 0 */
  int (*is_zero) (const polynomial_context_t* ctx, const coefficient_t* C);

  /** Returns true if the coefficient is 1 */
  int (*is_one) (const polynomial_context_t* ctx, const coefficient_t* C);

  /** Returns the sign of the coefficient in the model */
  int (*sgn) (const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m);

  /**
   * Returns the approximation of value of the coefficient (is either zero or
   * does not contain zero.
   */
  void (*value_approx) (const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m, interval_t* value);

  /** Returns the sign of the leading coefficient */
  int (*lc_sgn) (const polynomial_context_t* ctx, const coefficient_t* C);

  /** Returns treu if the coefficient is properly ordered according to ctx */
  int (*in_order) (const polynomial_context_t* ctx, const coefficient_t* C);

  /**
   * Compare the two coefficients in the ring. Not necessarily +/- 1, could be
   * any integer, only the sign matters.
   */
  int (*cmp) (const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

  /**
   * Compare the type of two coefficients. Type if either a constant or a
   * polynomial, where comparison is done only by the lead variable comparison
   * and constants are smaller than other polynomials. This is not a total
   * ORDER.
   */
  int (*cmp_type) (const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

  /** Returns true if C1 divides C2. */
  int (*divides) (const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

  /** Prints the coefficient to the given stream. */
  int (*print) (const polynomial_context_t* ctx, const coefficient_t* C, FILE* out);

  /** Returns the string representation of the coefficient. */
  char* (*to_string) (const polynomial_context_t* ctx, const coefficient_t* C);

  /** Order the polynomial according to the given order */
  void (*order) (const polynomial_context_t* ctx, coefficient_t* C);

  /** Compute S = C1 + C2. */
  void (*add) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute S = C1 - C2 in the given ring. */
  void (*sub) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute N = -C. */
  void (*neg) (const polynomial_context_t* ctx, coefficient_t* N, const coefficient_t* C);

  /** Compute P = C1 * C2. */
  void (*mul) (const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute P = C * a. */
  void (*mul_int) (const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, long a);

  /** Multiplication with x^n (x should be equal or biggier then top variable of C)  */
  void (*shl) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, variable_t x, unsigned n);

  /** Division with x^n (C should have degree at least n)  */
  void (*shr) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, unsigned n);

  /** Compute P = C^n. */
  void (*pow) (const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, unsigned n);

  /** Compute S += C1*C2. */
  void (*add_mul) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute S -= C1*C2. */
  void (*sub_mul) (const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

  /**
   * Reduce the coefficient A in Z[y,x] using B in Z[y,x] so that
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
   *   deg(R) < deg(B) or deg(R) = 0.
   *
   * If the exact flag is on, we assume that all division is exact.
   */
  void (*reduce) (const polynomial_context_t* ctx, const coefficient_t* A, const coefficient_t* B, coefficient_t* P, coefficient_t* Q, coefficient_t* R, remaindering_type_t type);

  /** Compute C1 = D*C2, in the given ring (assumes that C2 divides C1). */
  void (*div) (const polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute R = D*C2 - C1, in the given ring (assumes division is exact). */
  void (*rem) (const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute the pseudo remainder */
  void (*prem) (const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute the sparse pseudo remainder */
  void (*sprem) (const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

  /** Compute C1 = D*C2 + R, in the given ring (assumes division is exact). */
  void (*divrem) (const polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

  /** Computes the derivative of the coefficient (in the main variable) */
  void (*derivative) (const polynomial_context_t* ctx, coefficient_t* C_d, const coefficient_t* C);

  /** Compute the greatest common divisor gcd(C1, C2). */
  void (*gcd) (const polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2);

  /**
   * Compute the primitive part pp(A) (divide with content).
   */
  void (*pp) (const polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C);

  /**
   * Compute the content cont(A) (gcd of coefficients). Content is computed so that
   * A/cont(A) has a positive leading coefficient.
   */
  void (*cont) (const polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C);

  /** Comput pp and cont at the same time */
  void (*pp_cont) (const polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C);

  /** Compute the least common multiple. */
  void (*lcm) (const polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2);

  /**
   * Compute the resultant of C1 and C2 over their (common) top variable.
   */
  void (*resultant) (const polynomial_context_t* ctx, coefficient_t* res, const coefficient_t* C1, const coefficient_t* C2);

  /**
   * Compute the principal subresultant coefficients. PSC should be a constructed
   * array of size deg(B) + 1 filled with 0.
   */
  void (*psc) (const polynomial_context_t* ctx, coefficient_t* psc, const coefficient_t* C1, const coefficient_t* C2);

  /** Set the power symbol for print-outs */
  void (*set_power_symbol) (const char* pow);

} coefficient_ops_t;

/** Implementation of the coefficient operations */
extern const coefficient_ops_t coefficient_ops;
