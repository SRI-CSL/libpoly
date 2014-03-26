/*
 * polynomial_coefficient.h
 *
 *  Created on: Feb 4, 2014
 *      Author: dejan
 */

#pragma once

#include <polynomial_context.h>

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

typedef struct polynomial_rec_struct polynomial_rec_t;
typedef struct coefficient_struct coefficient_t;

/** Recursive nodes in the tree representation of the polynomial */
struct polynomial_rec_struct {
  /** The used size of the coefficient array */
  size_t size;
  /** Capacity of the coefficient array */
  size_t capacity;
  /** The main variable */
  variable_t x;
  /** Coefficients */
  coefficient_t* coefficients;
};

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
struct coefficient_struct {
  coefficient_type_t type;
  coefficient_union_t value;
};

#define SIZE(C) ((C)->value.rec.size)
#define CAPACITY(C) ((C)->value.rec.capacity)
#define COEFF(C, i) ((C)->value.rec.coefficients + (i))
#define VAR(C) ((C)->value.rec.x)

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

/** Construct a zero coefficient */
void coefficient_construct(const polynomial_context_t* ctx, coefficient_t* C);

/** Construct a coefficient from integer */
void coefficient_construct_from_int(const polynomial_context_t* ctx, coefficient_t* C, long C_int);

/** Construct a coefficient from integer */
void coefficient_construct_from_integer(const polynomial_context_t* ctx, coefficient_t* C, const integer_t* C_integer);

/** Construct from a univariate polynomial */
void coefficient_construct_from_univariate(const polynomial_context_t* ctx, coefficient_t* C, const upolynomial_t* p, variable_t x);

/** Construct a simple monomial coefficient a*x^n */
void coefficient_construct_simple(const polynomial_context_t* ctx, coefficient_t* C, const integer_t* a, variable_t x, unsigned n);

/** Construct a copy of the given coefficient. */
void coefficient_construct_copy(const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

/** Destructs the coefficient. */
void coefficient_destruct(coefficient_t* C);

/** Swap two coefficients. */
void coefficient_swap(coefficient_t* C1, coefficient_t* C2);

/** Assign the coefficient a given coefficient. */
void coefficient_assign(const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

/** Assign the coefficient a give integer */
void coefficient_assign_int(const polynomial_context_t* ctx, coefficient_t* C, long x);

/** Returns the univariate version of this coefficient (must be univariate) */
upolynomial_t* coefficient_to_univariate(const polynomial_context_t* ctx, const coefficient_t* C);

/** Returns true if the coefficient is a constant */
int coefficient_is_constant(const coefficient_t* C);

/** The degree of the coefficient */
size_t coefficient_degree(const coefficient_t* C);

/** The top variable of the coefficient */
variable_t coefficient_top_variable(const coefficient_t* C);

/** Get the k-th coefficient of this coefficient */
const coefficient_t* coefficient_get_coefficient(const coefficient_t* C, size_t k);

/** Get the leading coefficient of this coefficient */
const coefficient_t* coefficient_lc(const coefficient_t* C);

/** Returns true if the coefficient is 0 */
int coefficient_is_zero(const polynomial_context_t* ctx, const coefficient_t* C);

/** Returns true if the coefficient is 1 */
int coefficient_is_one(const polynomial_context_t* ctx, const coefficient_t* C);

/** Returns the sign of the coefficient in the model */
int coefficient_sgn(const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m);

/**
 * Returns the approximation of value of the coefficient (is either zero or
 * does not contain zero.
 */
void coefficient_value_approx(const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m, interval_t* value);

/** Returns the sign of the leading coefficient */
int coefficient_lc_sgn(const polynomial_context_t* ctx, const coefficient_t* C);

/** Returns treu if the coefficient is properly ordered according to ctx */
int coefficient_in_order(const polynomial_context_t* ctx, const coefficient_t* C);

/**
 * Compare the two coefficients in the ring. Not necessarily +/- 1, could be
 * any integer, only the sign matters.
 */
int coefficient_cmp(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compare the type of two coefficients. Type if either a constant or a
 * polynomial, where comparison is done only by the lead variable comparison
 * and constants are smaller than other polynomials. This is not a total
 * ORDER.
 */
int coefficient_cmp_type(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/** Returns true if C1 divides C2. */
int coefficient_divides(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/** Prints the coefficient to the given stream. */
int coefficient_print(const polynomial_context_t* ctx, const coefficient_t* C, FILE* out);

/** Returns the string representation of the coefficient. */
char* coefficient_to_string(const polynomial_context_t* ctx, const coefficient_t* C);

/** Order the polynomial according to the given order */
void coefficient_order(const polynomial_context_t* ctx, coefficient_t* C);

/** Compute S = C1 + C2. */
void coefficient_add(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute S = C1 - C2 in the given ring. */
void coefficient_sub(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute N = -C. */
void coefficient_neg(const polynomial_context_t* ctx, coefficient_t* N, const coefficient_t* C);

/** Compute P = C1 * C2. */
void coefficient_mul(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C1, const coefficient_t* C2);

/** Compute P = C * a. */
void coefficient_mul_int(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, long a);

/** Multiplication with x^n (x should be equal or biggier then top variable of C)  */
void coefficient_shl(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, variable_t x, unsigned n);

/** Division with x^n (C should have degree at least n)  */
void coefficient_shr(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, unsigned n);

/** Compute P = C^n. */
void coefficient_pow(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, unsigned n);

/** Compute S += C1*C2. */
void coefficient_add_mul(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute S -= C1*C2. */
void coefficient_sub_mul(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

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
void coefficient_reduce(const polynomial_context_t* ctx, const coefficient_t* A, const coefficient_t* B, coefficient_t* P, coefficient_t* Q, coefficient_t* R, remaindering_type_t type);

/** Compute C1 = D*C2, in the given ring (assumes that C2 divides C1). */
void coefficient_div(const polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2);

/** Divide the degrees of the main variable of coefficient with the given number */
void coefficient_div_degrees(const polynomial_context_t* ctx, coefficient_t* C, size_t p);

/** Compute R = D*C2 - C1, in the given ring (assumes division is exact). */
void coefficient_rem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute the pseudo remainder */
void coefficient_prem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute the sparse pseudo remainder */
void coefficient_sprem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute C1 = D*C2 + R, in the given ring (assumes division is exact). */
void coefficient_divrem(const polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Computes the derivative of the coefficient (in the main variable) */
void coefficient_derivative(const polynomial_context_t* ctx, coefficient_t* C_d, const coefficient_t* C);

/** Compute the greatest common divisor gcd(C1, C2). */
void coefficient_gcd(const polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the primitive part pp(A(divide with content).
 */
void coefficient_pp(const polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C);

/**
 * Compute the content cont(A(gcd of coefficients). Content is computed so that
 * A/cont(A) has a positive leading coefficient.
 */
void coefficient_cont(const polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C);

/** Comput pp and cont at the same time */
void coefficient_pp_cont(const polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C);

/** Compute the least common multiple. */
void coefficient_lcm(const polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the resultant of C1 and C2 over their (common) top variable.
 */
void coefficient_resultant(const polynomial_context_t* ctx, coefficient_t* res, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the principal subresultant coefficients. PSC should be a constructed
 * array of size deg(B) + 1 filled with 0.
 */
void coefficient_psc(const polynomial_context_t* ctx, coefficient_t* psc, const coefficient_t* C1, const coefficient_t* C2);

/** Set the power symbol for print-outs */
void coefficient_set_power_symbol(const char* pow);

