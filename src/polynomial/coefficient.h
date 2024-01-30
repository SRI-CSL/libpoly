/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <polynomial_context.h>
#include <monomial.h>
#include <assignment.h>

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
  lp_variable_t x;
  /** Coefficients */
  coefficient_t* coefficients;
};

/**
 * Value of each coefficient is either a base value or a polynomial.
 */
typedef union {
  lp_integer_t num;
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
void coefficient_construct(const lp_polynomial_context_t* ctx, coefficient_t* C);

/** Construct a coefficient from integer */
void coefficient_construct_from_int(const lp_polynomial_context_t* ctx, coefficient_t* C, long C_int);

/** Construct a coefficient from integer */
void coefficient_construct_from_integer(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* C_integer);

/** Construct from a univariate polynomial */
void coefficient_construct_from_univariate(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_upolynomial_t* p, lp_variable_t x);

/** Construct a simple monomial coefficient a*x^n */
void coefficient_construct_simple(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* a, lp_variable_t x, unsigned n);

/** Construct a simple monomial coefficient a*x^n */
void coefficient_construct_simple_int(const lp_polynomial_context_t* ctx, coefficient_t* C, long a, lp_variable_t x, unsigned n);

/** Construcut a simple linear polynomial a*x + b */
void coefficient_construct_linear(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* a, const lp_integer_t* b, lp_variable_t x);

/** Construct a copy of the given coefficient. */
void coefficient_construct_copy(const lp_polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

/** Destructs the coefficient. */
void coefficient_destruct(coefficient_t* C);

/** Swap two coefficients. */
void coefficient_swap(coefficient_t* C1, coefficient_t* C2);

/** Assign the coefficient a given coefficient. */
void coefficient_assign(const lp_polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from);

/** Assign the coefficient a give integer */
void coefficient_assign_int(const lp_polynomial_context_t* ctx, coefficient_t* C, long x);

/** Assign the coefficient a give integer */
void coefficient_assign_integer(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* x);

/** Check if the coefficient is univariate */
int coefficient_is_univariate(const coefficient_t* C);

/** Check if the coefficient is linear */
int coefficient_is_linear(const coefficient_t* C);

/** Returns true if the coefficient conforms to internal representation (mainly debug purposes) */
int coefficient_is_normalized(const lp_polynomial_context_t* ctx, coefficient_t* C);

/** Returns the univariate version of this coefficient (must be univariate) */
lp_upolynomial_t* coefficient_to_univariate(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** Returns the univariate version of this coefficient with regard to m. m must evaluate all (sub-)coefficients of C to integer. */
lp_upolynomial_t* coefficient_to_univariate_m(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m);

/** Returns true if the coefficient is a constant */
int coefficient_is_constant(const coefficient_t* C);

/** Returns the constant of the polynomial (the actual deep constant) */
const lp_integer_t* coefficient_get_constant(const coefficient_t* C);

/** Returns true if the coefficient is a monomial */
int coefficient_is_monomial(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** In case C is a monomial, returns it. */
void coefficient_to_monomial(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_monomial_t *out);

/** The degree of the coefficient */
size_t coefficient_degree(const coefficient_t* C);

/** The degree of the coefficient in the model */
size_t coefficient_degree_m(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M);

/** Degree over variable x (only call if x is bigger than rest) */
size_t coefficient_degree_safe(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_variable_t x);

/** The top variable of the coefficient */
lp_variable_t coefficient_top_variable(const coefficient_t* C);

/** Get the k-th coefficient of this coefficient */
const coefficient_t* coefficient_get_coefficient(const coefficient_t* C, size_t k);

/** Get the k-th coefficient of this coefficient */
const coefficient_t* coefficient_get_coefficient_safe(const lp_polynomial_context_t* ctx, const coefficient_t* C, size_t k, lp_variable_t x);

/** Get the leading coefficient of this coefficient */
const coefficient_t* coefficient_lc(const coefficient_t* C);

/** Get the leading coefficient of this coefficient (only when x top or bigger) */
const coefficient_t* coefficient_lc_safe(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_variable_t x);

/** Get the model-based leading coefficient */
const coefficient_t* coefficient_lc_m(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M);

/** Get the reductum of the coefficient */
void coefficient_reductum(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C);

/** Get the model-based reductum of the coefficient */
void coefficient_reductum_m(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C, const lp_assignment_t* m, lp_polynomial_vector_t* assumptions);

/** Returns true if the coefficient is 0 */
int coefficient_is_zero(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** Returns true if the coefficient is 1 */
int coefficient_is_one(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** Returns true if the coefficient is -1 */
int coefficient_is_minus_one(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** Returns true if all variables of C are assigned */
int coefficient_is_assigned(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m);

/** Returns the sign of the coefficient in the model */
int coefficient_sgn(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m);

/** Returns the interval approximation of the value of the polynomial */
void coefficient_interval_value(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_interval_assignment_t* m, lp_interval_t* result);


/**
 * Returns the approximation of value of the coefficient (is either zero or
 * does not contain zero.
 */
void coefficient_value_approx(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m, lp_rational_interval_t* value);

/** Evaluates the coefficient in the model. */
lp_value_t* coefficient_evaluate(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M);

/** Returns the sign of the leading coefficient */
int coefficient_lc_sgn(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/** Returns the integer value of the leading coefficient */
void coefficient_lc_constant(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_integer_t* out);

/** Returns true if the coefficient is properly ordered according to ctx */
int coefficient_in_order(const lp_polynomial_context_t* ctx, const coefficient_t* C);

/**
 * Compare the two coefficients in the ring. Not necessarily +/- 1, could be
 * any integer, only the sign matters.
 */
int coefficient_cmp(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compare the type of two coefficients. Type if either a constant or a
 * polynomial, where comparison is done only by the lead variable comparison
 * and constants are smaller than other polynomials. This is not a total
 * ORDER.
 */
int coefficient_cmp_type(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/** Returns true if C1 divides C2. */
int coefficient_divides(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2);

/** Order the polynomial according to the given order */
void coefficient_order(const lp_polynomial_context_t* ctx, coefficient_t* C);

/** Compute S = C1 + C2. */
void coefficient_add(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute S = C1 - C2 in the given ring. */
void coefficient_sub(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute N = -C. */
void coefficient_neg(const lp_polynomial_context_t* ctx, coefficient_t* N, const coefficient_t* C);

/** Compute P = C1 * C2. */
void coefficient_mul(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C1, const coefficient_t* C2);

/** Compute P = C * a. */
void coefficient_mul_integer(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, const lp_integer_t* a);

/** Compute P = C * a. */
void coefficient_mul_int(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, long a);

/** Multiplication with x^n (x should be equal or bigger then top variable of C)  */
void coefficient_shl(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, lp_variable_t x, unsigned n);

/** Division with x^n (C should have degree at least n)  */
void coefficient_shr(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, lp_variable_t x, unsigned n);

/** Compute P = C^n. */
void coefficient_pow(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, unsigned n);

/** Compute S += C1*C2. */
void coefficient_add_mul(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/** Compute S -= C1*C2. */
void coefficient_sub_mul(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Reduce the coefficient A in Z[x,y] using B in Z[x,y] so that
 *
 *   P*A = Q*B + R
 *
 * and
 *
 *   P in Z[x]
 *   Q, R in Z[x,y]
 *
 * with
 *
 *   deg(R) < deg(B) or deg(R) = 0.
 *
 * If the exact flag is on, we assume that all division is exact.
 */
void coefficient_reduce(const lp_polynomial_context_t* ctx, const coefficient_t* A, const coefficient_t* B, coefficient_t* P, coefficient_t* Q, coefficient_t* R, remaindering_type_t type);

/** Compute C1 = D*C2, in the given ring (assumes that C2 divides C1). */
void coefficient_div(const lp_polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2);

/** Compute C = C/A */
void coefficient_div_constant(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* A);

/** Divide the degrees of the main variable of coefficient with the given number */
void coefficient_div_degrees(const lp_polynomial_context_t* ctx, coefficient_t* C, size_t p);

/** Divides the degree of exponents wrt to the prime field order */
void coefficient_reduce_Zp(const lp_polynomial_context_t* ctx, coefficient_t* C);

/** Compute R = D*C2 - C1, in the given ring (assumes division is exact). */
void coefficient_rem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute the pseudo remainder */
void coefficient_prem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute a*C1 = D*C2 + R, in the given ring (using pseudo division). */
void coefficient_pdivrem(const lp_polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute the sparse pseudo remainder */
void coefficient_sprem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute a*C1 = D*C2 + R, in the given ring (using sparse pseudo division). */
void coefficient_spdivrem(const lp_polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Compute C1 = D*C2 + R, in the given ring (assumes division is exact). */
void coefficient_divrem(const lp_polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2);

/** Computes the derivative of the coefficient (in the main variable) */
void coefficient_derivative(const lp_polynomial_context_t* ctx, coefficient_t* C_d, const coefficient_t* C);

/**
 * Compute the resultant of C1 and C2 over their (common) top variable.
 */
void coefficient_resultant(const lp_polynomial_context_t* ctx, coefficient_t* res, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the subresultant regular sub-chain.
 * SRS should be a constructed array of size deg(B) + 1 filled with 0.
 */
void coefficient_srs(const lp_polynomial_context_t* ctx, coefficient_t* srs, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the principal subresultant coefficients. These are the leading coefficients of the SRS.
 * PSC should be a constructed array of size deg(B) + 1 filled with 0.
 */
void coefficient_psc(const lp_polynomial_context_t* ctx, coefficient_t* psc, const coefficient_t* C1, const coefficient_t* C2);

/** Function type called on coefficient traversal */
typedef void (*traverse_f) (const lp_polynomial_context_t* ctx, lp_monomial_t* p, void* data);

/**
 * Run on the coefficient to traverse all monomials. The traverse_f function will be called on he monomial with
 * the associated data.
 */
void coefficient_traverse(const lp_polynomial_context_t* ctx, const coefficient_t* C, traverse_f f, lp_monomial_t* m, void* data);

/**
 * Method called to add a monomial to C. The monomial should be ordered in the
 * same order as C, top variable at the m[0].
 */
void coefficient_add_ordered_monomial(const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* C_void);

/**
 * Method called to add a monomial to C.
 */
void coefficient_order_and_add_monomial(const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* C_void);

/**
 * Add the monomial to the coefficient.
 */
void coefficient_add_monomial(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_monomial_t* m);

/**
 * Substitute and evaluate any integer value. M must assign to integer only.
 * Can handle ctx->K != lp_Z.
 */
void coefficient_evaluate_integer(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M, lp_integer_t *out);

/**
 * Substitute and evaluate any rational values. The result C_out = multiplier*M(C) with the multiplier
 * being > 0.
 */
void coefficient_evaluate_rationals(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M, coefficient_t* C_out, lp_integer_t* multiplier);

/**
 * Get the variables of the coefficient.
 */
void coefficient_get_variables(const coefficient_t* C, lp_variable_list_t* vars);

/**
 * Isolate the roots (multivariate with model).
 */
void coefficient_roots_isolate(const lp_polynomial_context_t* ctx, const coefficient_t* A, const lp_assignment_t* M, lp_value_t* roots, size_t* roots_size);

/**
 * Isolate the roots (univaraite no model).
 */
void coefficient_roots_isolate_univariate(const lp_polynomial_context_t* ctx, const coefficient_t* A, lp_value_t* roots, size_t* roots_size);

/**
 * Get the hash of the polynomial. The hash is simple to ensure we keep the hash across different
 * orders.
 */
size_t coefficient_hash(const lp_polynomial_context_t* ctx, const coefficient_t* A);
