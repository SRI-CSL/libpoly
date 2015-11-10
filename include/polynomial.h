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

#include "poly.h"

#include "variable_order.h"
#include "integer.h"
#include "polynomial_context.h"
#include "assignment.h"
#include "monomial.h"
#include "sign_condition.h"

/**
 * Polynomials incorporate the context and coefficient data. The also carry
 * flags so as not to re-do any computation.
 *
 * If the external flag is on during construction, the polynomial
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

/** Construct a zero polynomial (does not attach the context) */
void lp_polynomial_construct(lp_polynomial_t* A, const lp_polynomial_context_t* ctx);

/** Construct a simple polynomial c*x^n */
void lp_polynomial_construct_simple(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const lp_integer_t* c, lp_variable_t x, unsigned n);

/** Construct a copy of the given polynomial (does not attach the context). */
void lp_polynomial_construct_copy(lp_polynomial_t* A, const lp_polynomial_t* from);

/** Destruct the polynomial. */
void lp_polynomial_destruct(lp_polynomial_t* A);

/** Delete the polynomial (destrcuct and free) */
void lp_polynomial_delete(lp_polynomial_t* A);

/** Allocate a new polynomial (unconstructed) */
lp_polynomial_t* lp_polynomial_alloc(void);

/** Allocate and construct a new polynomial */
lp_polynomial_t* lp_polynomial_new(const lp_polynomial_context_t* ctx);

/** Allocate and construct a copy of the given polynomial */
lp_polynomial_t* lp_polynomial_new_copy(const lp_polynomial_t* A);

/** Make the polynomial as external */
void lp_polynomial_set_external(lp_polynomial_t* A);

/** Swap two polynomials. */
void lp_polynomial_swap(lp_polynomial_t* A1, lp_polynomial_t* A2);

/** Assign the polynomial a given polynomial. */
void lp_polynomial_assign(lp_polynomial_t* A, const lp_polynomial_t* from);

/** Returns the context of the polynomial */
const lp_polynomial_context_t* lp_polynomial_context(const lp_polynomial_t* A);

/** Returns the degree of the polynomial (in it's top variable) */
size_t lp_polynomial_degree(const lp_polynomial_t* A);

/** Returns the top variable of the polynomial */
lp_variable_t lp_polynomial_top_variable(const lp_polynomial_t* A);

/** Returns 1 if the leading coefficient is constant */
int lp_polynomial_lc_is_constant(const lp_polynomial_t* A);

/** In case lc is constant, this returns the sign */
int lp_polynomial_lc_sgn(const lp_polynomial_t* A);

/** Get the context of the given polynomial */
const lp_polynomial_context_t* lp_polynomial_get_context(const lp_polynomial_t* A);

/** Returns all the variables (it will not clear the output list vars) */
void lp_polynomial_get_variables(const lp_polynomial_t* A, lp_variable_list_t* vars);

/** Puts the k-th coefficient of A into C */
void lp_polynomial_get_coefficient(lp_polynomial_t* C, const lp_polynomial_t* A, size_t k);

/** Get the reductum of the polynomial (the polynomial without the leading coefficient) */
void lp_polynomial_reductum(lp_polynomial_t* R, const lp_polynomial_t* A);

/** Get the model-based reductum of the polynomial (the polynomial without the leading coefficient) */
void lp_polynomial_reductum_m(lp_polynomial_t* R, const lp_polynomial_t* A, const lp_assignment_t* m);

/** Returns true if the polynomial is a constant */
int lp_polynomial_is_constant(const lp_polynomial_t* A);

/** Returns true if the polynomial is 0 */
int lp_polynomial_is_zero(const lp_polynomial_t* A);

/** Returns true if the polynomial is univariate */
int lp_polynomial_is_univariate(const lp_polynomial_t* A);

/**
 * Returns true if polynomial is univariate, i.e. has one variable unassigned,
 * and the unassigned variable is the top variable.
 */
int lp_polynomial_is_univariate_m(const lp_polynomial_t* A, const lp_assignment_t* m);

/** Returns the univariate polynomial (if univariate, or 0 otherwise) */
lp_upolynomial_t* lp_polynomial_to_univariate(const lp_polynomial_t* A);

/** returns the sign of the polynomial in the model (-1, 0, +1), or -2 if not all variables assigned */
int lp_polynomial_sgn(const lp_polynomial_t* A, const lp_assignment_t* m);

/** returns the sign of the polynomial in the model */
lp_value_t* lp_polynomial_evaluate(const lp_polynomial_t* A, const lp_assignment_t* m);


/**
 * Compare the two polynomials in the ring. Not necessarily +/- 1, could be
 * any integer, only the sign matters.
 */
int lp_polynomial_cmp(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/**
 * Compare the type of two polynomials. Type if either a constant or a
 * polynomial, where comparison is done only by the lead variable comparison
 * and constants are smaller than other polynomials. This is not a total
 * ORDER.
 */
int lp_polynomial_cmp_type(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Returns the hash of the polynomial */
size_t lp_polynomial_hash(const lp_polynomial_t* A);

/** Compares two polynomials for equality. */
int lp_polynomial_eq(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Returns true if A1 divides A2. */
int lp_polynomial_divides(const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Prints the polynomial to the given stream. */
int lp_polynomial_print(const lp_polynomial_t* A, FILE* out);

#if HAVE_OPEN_MEMSTREAM
/** Returns the string representation of the polynomial. */
char* lp_polynomial_to_string(const lp_polynomial_t* A);
#endif

/** Compute S = A1 + A2. */
void lp_polynomial_add(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute S += M. */
void lp_polynomial_add_monomial(lp_polynomial_t* S, const lp_monomial_t* M);

/** Compute S = A1 - A2 in the given ring. */
void lp_polynomial_sub(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute N = -C. */
void lp_polynomial_neg(lp_polynomial_t* N, const lp_polynomial_t* A);

/** Compute P = A1 * A2. */
void lp_polynomial_mul(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute P = A1 * C. */
void lp_polynomial_mul_integer(lp_polynomial_t* P, const lp_polynomial_t* A1, const lp_integer_t* C);

/**
 * Multiplication with x^n.
 */
void lp_polynomial_shl(lp_polynomial_t* S, const lp_polynomial_t* A, unsigned n);

/** Compute P = C^n. */
void lp_polynomial_pow(lp_polynomial_t* P, const lp_polynomial_t* A, unsigned n);

/** Compute S += A1*A2. */
void lp_polynomial_add_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute S -= A1*A2. */
void lp_polynomial_sub_mul(lp_polynomial_t* S, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

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
void lp_polynomial_reduce(const lp_polynomial_t* A, const lp_polynomial_t* B,
    lp_polynomial_t* P, lp_polynomial_t* Q, lp_polynomial_t* R);

/** Compute A1 = D*A2 (assumes that A2 divides A1). */
void lp_polynomial_div(lp_polynomial_t* D, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute A1 = D*A2 + R (assumes that exact division). */
void lp_polynomial_rem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute a*A1 = D*A2 + R (pseudo remainder). */
void lp_polynomial_prem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute a*A1 = D*A2 + R (sparse pseudo remainder). */
void lp_polynomial_sprem(lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute A1 = D*A2 + R (assumes that exact division). */
void lp_polynomial_divrem(lp_polynomial_t* D, lp_polynomial_t* R, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute A_d = A' (in the top variable). */
void lp_polynomial_derivative(lp_polynomial_t* A_d, const lp_polynomial_t* A);

/** Compute the greatest common divisor gcd(A1, A2). */
void lp_polynomial_gcd(lp_polynomial_t* gcd, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/** Compute the least common multiple lcm(A1, A2). */
void lp_polynomial_lcm(lp_polynomial_t* lcm, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/**
 * Compute the resultant of A1 and A2 in their top variable. Both A1 and A2
 * must be (non-trivial) polynomials over the same variable.
 */
void lp_polynomial_resultant(lp_polynomial_t* res, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/**
 * Compute the principal subresultant coefficients (psc) of A1 and A1. Bot A1
 * and A2 must be (non-trivial) polynomials over the same variable. If
 *  deg(A1) = m, deg(A2) = n, and output will be of size min(m, n).
 */
void lp_polynomial_psc(lp_polynomial_t** psc, const lp_polynomial_t* A1, const lp_polynomial_t* A2);

/**
 * Get the square-free factorization of the given polynomial. It will allocate
 * the given arrays and return the polynomials and their multiplicities in
 * them.
 */
void lp_polynomial_factor_square_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size);

/**
 * Get the content-free factorization of the given polynomial. It will allocate
 * the given arrays and return the polynomials and their multiplicities in
 * them.
 */
void lp_polynomial_factor_content_free(const lp_polynomial_t* A, lp_polynomial_t*** factors, size_t** multiplicities, size_t* size);

/**
 * Given a polynomial p(x1, ..., xn, y) with y being the top variable, and an
 * assignment M that assigns x1, ..., xn, the function returns the roots of p
 * at M. The output array should be big enough to fit the root, i.e. at least
 * the degree of y. The array should be unconstructed, but you should destruct it
 * once done.
 */
void lp_polynomial_roots_isolate(const lp_polynomial_t* A, const lp_assignment_t* M, lp_value_t* roots, size_t* roots_size);

/**
 * Given a polynomial A(x1, ..., xn, y) with y being the top variable, a sign
 * condition, and an assignment M that assigns x1, ..., xn, the function returns
 * a subset of R where
 *
 *   sgn(A(M(x1), ..., M(xn), y)) = sgn_condition .
 *
 * If negated is true, the constraint is considered negated.
 */
lp_feasibility_set_t* lp_polynomial_constraint_get_feasible_set(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, int negated, const lp_assignment_t* M);

/**
 * Given a polynomial constraint, as above, evaluate its truth value.
 */
int lp_polynomial_constraint_evaluate(const lp_polynomial_t* A, lp_sign_condition_t sgn_condition, const lp_assignment_t* M);

/**
 * Given a polynomial A(x1, ..., xn, y) with y being the top variable, a root index,
 * a sign condition, and an assignment M that assigns x1, ..., xn, the function
 * returns a subset or R where
 *
 * root_index < root_cound(A, M) && sgn(x - root(k, A(M(x1, ..., M(xn), y)) == sgn_condition
 *
 * If negated is true, the constraint is considered negated.
 */
lp_feasibility_set_t* lp_polynomial_root_constraint_get_feasible_set(const lp_polynomial_t* A, size_t root_index, lp_sign_condition_t sgn_condition, int negated, const lp_assignment_t* M);

/**
 * Given a root constraint as above, evaluate its truth value.
 */
int lp_polynomial_root_constraint_evaluate(const lp_polynomial_t* A, size_t root_index, lp_sign_condition_t sgn_condition, const lp_assignment_t* M);

/**
 * Function type called on polynomial traversal. It will be called on all monomials.
 */
typedef void (*lp_polynomial_traverse_f) (const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* data);

/**
 * Run on the polynomial to traverse all monomials. The traverse_f function will be called on the monomial with
 * the associated data.
 */
void lp_polynomial_traverse(const lp_polynomial_t* A, lp_polynomial_traverse_f f, void* data);

/** Check the intergrity of the polynomial */
int lp_polynomial_check_integrity(const lp_polynomial_t* A);
