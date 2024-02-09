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

#include "integer.h"
#include "rational.h"
#include "interval.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Construct the polynomial given its coefficients. Coefficients should be
 * indexed by degree and they will be normalized according to the given ring.
 */
lp_upolynomial_t* lp_upolynomial_construct(const lp_int_ring_t* K, size_t degree, const lp_integer_t* coefficients);

/**
 * Construct the polynomial c*x^d.
 */
lp_upolynomial_t* lp_upolynomial_construct_power(const lp_int_ring_t* K, size_t degree, long c);

/**
 * Construct the polynomial given its coefficients. Coefficients should be
 * indexed by degree and they will be normalize according to the given ring.
 */
lp_upolynomial_t* lp_upolynomial_construct_from_int(const lp_int_ring_t* K, size_t degree, const int* coefficients);

/**
 * Construct the polynomial given its coefficients. Coefficients should be
 * indexed by degree and they will be normalize according to the given ring.
 */
lp_upolynomial_t* lp_upolynomial_construct_from_long(const lp_int_ring_t* K, size_t degree, const long* coefficients);

/**
 * Construct a copy of the polynomial.
 */
lp_upolynomial_t* lp_upolynomial_construct_copy(const lp_upolynomial_t* p);

/**
 * Construct a copy of the polynomial, but change the ring.
 */
lp_upolynomial_t* lp_upolynomial_construct_copy_K(const lp_int_ring_t* K, const lp_upolynomial_t* p);

/**
 * Frees the polynomial data and detaches the ring. */
void lp_upolynomial_delete(lp_upolynomial_t* p);

/**
 * Returns the degree of the polynomial. Note that the degree of the constant
 * 0 is 0.
 */
size_t lp_upolynomial_degree(const lp_upolynomial_t* p);

/**
 * Returns the field of the polynomial (unatached).
 */
const lp_int_ring_t* lp_upolynomial_ring(const lp_upolynomial_t* p);

/**
 * Sets the ring to given ring (has to be "larger" than existing).
 */
void lp_upolynomial_set_ring(lp_upolynomial_t* p, const lp_int_ring_t* K);

/**
 * Returns the lead coefficient of the given polynomial.
 */
const lp_integer_t* lp_upolynomial_lead_coeff(const lp_upolynomial_t* p);

/**
 * Returns the constant term of the given polynomial, or 0 if no constant factor.
 */
const lp_integer_t* lp_upolynomial_const_term(const lp_upolynomial_t* p);

/**
 * Unpack the polynomial into a dense representation. The out vector is
 * assumed to be large enough, and filled with 0 (only non-zero coefficients
 * will be copied into out).
 */
void lp_upolynomial_unpack(const lp_upolynomial_t* p, lp_integer_t* out);

/**
 * Print the polynomial to the output stream.
 */
int lp_upolynomial_print(const lp_upolynomial_t* p, FILE* out);

/**
 * Get a string representation of the polynomial (you own the memory).
 */
char* lp_upolynomial_to_string(const lp_upolynomial_t* p);

/**
 * Converts the upolynomial to a polynomial in variable var.
 * ctx must have the same ring than upolynomial.
 */
lp_polynomial_t* lp_upolynomial_to_polynomial(const lp_upolynomial_t* p, const lp_polynomial_context_t* ctx, lp_variable_t var);

/**
 * Returns true if this is a zero polynomial
 */
int lp_upolynomial_is_zero(const lp_upolynomial_t* p);

/**
 * Returns true if this is the polynomial 1.
 */
int lp_upolynomial_is_one(const lp_upolynomial_t* p);

/**
 * Returns true if the polynomial is monic.
 */
int lp_upolynomial_is_monic(const lp_upolynomial_t* p);

/**
 * Returns true if the polynomial is primitive, i.e. gcd of all coefficients
 * is 1 and the leading coefficient is positive.
 */
int lp_upolynomial_is_primitive(const lp_upolynomial_t* p);

/**
 * Evaluates the polynomial.
 */
void lp_upolynomial_evaluate_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x, lp_integer_t* value);

/**
 * Evaluates the polynomial. Only makes sense for polynomials in Z[x].
 */
void lp_upolynomial_evaluate_at_rational(const lp_upolynomial_t* p, const lp_rational_t* x, lp_rational_t* value);

/**
 * Evaluates the polynomial. Only makes sense for polynomials in Z[x].
 */
void lp_upolynomial_evaluate_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x, lp_dyadic_rational_t* value);

/**
 * Get the sign of the polynomial in the given integer point.
 */
int lp_upolynomial_sgn_at_integer(const lp_upolynomial_t* p, const lp_integer_t* x);

/**
 * Get the sign of the polynomial in the given rational point. Only makes sense
 * for polynomials in Z[x].
 */
int lp_upolynomial_sgn_at_rational(const lp_upolynomial_t* p, const lp_rational_t* x);

/**
 * Get the sign of the polynomial in the given rational point. Only makes sense
 * for polynomials in Z[x].
 */
int lp_upolynomial_sgn_at_dyadic_rational(const lp_upolynomial_t* p, const lp_dyadic_rational_t* x);

/**
 * Compares two polynomials (lexicographic from highest coefficient) and
 * returns -1 if p < q, 0 if p == q, and 1 if p > q.
 */
int lp_upolynomial_cmp(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Returns the polynomial f(-x).
 */
lp_upolynomial_t* lp_upolynomial_subst_x_neg(const lp_upolynomial_t* f);

/**
 * Replace the polynomial with f(x**n).
 * The function suppose that the maximum degree * n doesn't overflow.
 */
void lp_upolynomial_subst_x_pow_in_place(lp_upolynomial_t* f, size_t n);

/**
 * Returns the polynomial -f.
 */
lp_upolynomial_t* lp_upolynomial_neg(const lp_upolynomial_t* p);

/**
 * Negates p in place.
 */
void lp_upolynomial_neg_in_place(lp_upolynomial_t* p);

/**
 * Makes the polynomial monic in place.
 */
lp_upolynomial_t* lp_upolynomial_make_monic(const lp_upolynomial_t* p);

/**
 * Makes the polynomial monic in place.
 */
void lp_upolynomial_make_monic_in_place(lp_upolynomial_t* p);

/**
 * Add two polynomials (all operations in the same ring).
 */
lp_upolynomial_t* lp_upolynomial_add(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Subtract two polynomials (all operations in the same ring).
 */
lp_upolynomial_t* lp_upolynomial_sub(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Multiply the two polynomials (all operations in the same ring).
 */
lp_upolynomial_t* lp_upolynomial_mul(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Multiply the two polynomials (all operations in the ring of p).
 */
lp_upolynomial_t* lp_upolynomial_mul_c(const lp_upolynomial_t* p, const lp_integer_t* c);

/**
 * Power of a polynomial (all operations in the ring of p).
 */
lp_upolynomial_t* lp_upolynomial_pow(const lp_upolynomial_t* p, long pow);

/**
 * Returns the derivative of the given polynomial. Note that deg(p') can be
 * less than deg(p) - 1 in some rings. For example in Z_4 (2x^2)' = 0.
 */
lp_upolynomial_t* lp_upolynomial_derivative(const lp_upolynomial_t* p);

/***
 * Returns true if p divides q.
 */
int lp_upolynomial_divides(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Return a polynomial with all the monomial degrees divided by the given
 * positive number a. All degrees must be divisible by a.
 */
lp_upolynomial_t* lp_upolynomial_div_degrees(const lp_upolynomial_t* p, size_t a);

/**
 * Returns the exact division of two polynomials. This assumes that p and q
 * are in the same ring and all needed coefficient inverses can be computed.
 */
lp_upolynomial_t* lp_upolynomial_div_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Returns the exact division of the polynomial with a constant c. This
 * assumes that all coefficients of p are divisible by c.
 */
lp_upolynomial_t* lp_upolynomial_div_exact_c(const lp_upolynomial_t* p, const lp_integer_t* c);

/**
 * Returns the exact remainder of two polynomials. This assumes that p and q
 * are in the same ring all needed coefficient inverses can be computed.
 */
lp_upolynomial_t* lp_upolynomial_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Returns the exact division of two polynomials. This assumes that p and q
 * are in the same ring and all needed coefficient inverses can be computed.
 * The output in div and rem will be newly allocated.
 */
void lp_upolynomial_div_rem_exact(const lp_upolynomial_t* p, const lp_upolynomial_t* q,
    lp_upolynomial_t** div, lp_upolynomial_t** rem);

/**
 * Pseudo-division of polynomials, div and rem such that
 *
 *   lc(q)^(p_deg - q_deg + 1) p = div*q + rem
 *
 * This assumes that deg(p) >= deg(q).
 *
 * Note: all computation is done in ring of p and q, but the algorithm doesn't
 * take advantage of existence of possible inverses -- algorithm proceeds as
 * if done in Z with individual operations performed in the ring.
 */
void lp_upolynomial_div_pseudo(lp_upolynomial_t** div, lp_upolynomial_t** rem, const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Compute the content of the polynomial. Content of the polynomial is the
 * gcd of the coefficients, of the same sign as the leading coefficient.
 * NOTE: The gcd is computed in Z so p must be in Z
 */
void lp_upolynomial_content_Z(const lp_upolynomial_t* p, lp_integer_t* content);

/**
 * Make the polynomial primitive. The polynomial is primitive if the content
 * is 1.
 */
void lp_upolynomial_make_primitive_Z(lp_upolynomial_t* p);

/**
 * Get the primitive part of the polynomial, the primitive part is p/content
 * and always has leading coefficient > 0.
 */
lp_upolynomial_t* lp_upolynomial_primitive_part_Z(const lp_upolynomial_t* p);

/**
 * Computes the polynomial greatest common divisor of p and q. The rings of p
 * and q are assumed to be the same. The lc(gcd) > 0, and if the ring of p
 * and q is a field, it will be normalized to monic -- lc(gcd) == 1.
 */
lp_upolynomial_t* lp_upolynomial_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q);

/**
 * Computes the extended gcd of p and q, i.e.
 *
 *   u*p + v*q = gcd(p, q).
 *
 * The coefficient rings of p and q are assumed to be the same. In addition
 * they are assumed to be prime fields. The function ensures a monic
 * GCD -- lc(gcd) == 1.
 *
 * All resulting polynomials (gcd, *u, *v) will be freshly allocated.
 */
lp_upolynomial_t* lp_upolynomial_extended_gcd(const lp_upolynomial_t* p, const lp_upolynomial_t* q, lp_upolynomial_t** u, lp_upolynomial_t** v);

/**
 * Given p, q, and r solve the equation
 *
 *  u*p+v*q = r
 *
 * for u and v. Assumes that gcd(p, q) divides r. Result such that
 *
 *   deg(u) < deg(q), deg(v) < deg(p)
 *
 * The coefficient rings of p, q, and r, are assumed to be the same. In
 * addition they are assumed to be a prime field.
 */
void lp_upolynomial_solve_bezout(const lp_upolynomial_t* p, const lp_upolynomial_t* q, const lp_upolynomial_t* r,
    lp_upolynomial_t** u, lp_upolynomial_t** v);

/**
 * Returns the factorization of the given polynomial in its ring. (DO NOT USE)
 */
lp_upolynomial_factors_t* lp_upolynomial_factor(const lp_upolynomial_t* p);

/**
 * Returns the square-free factorization of the given polynomial in its ring.
 * In a square-free factorization each factor is square-free. Individual
 * factors are also mutually prime, i.e. gcd(f_i, f_j) = 1 for i != j. If x
 * appears as a factor, it is always separate in the result.
 */
lp_upolynomial_factors_t* lp_upolynomial_factor_square_free(const lp_upolynomial_t* p);

/**
 * Return the Sturm sequence of the given polynomial. The arrays S will be
 * allocated, and the user should de-allocate it. The size parameter will be
 * updated with the size of the array.
 */
void lp_upolynomial_sturm_sequence(const lp_upolynomial_t* f, lp_upolynomial_t*** S, size_t* size);

/**
 * Counts the number of real roots in the given interval. If the interval is
 * 0, it counts through (-inf, inf).
 */
int lp_upolynomial_roots_count(const lp_upolynomial_t* p, const lp_rational_interval_t* ab);

/**
 * Isolate the distinct real roots of the given polynomial. Roots should be
 * an (unconstructed) array of algebraic numbers with size at least
 * lp_upolynomial_roots_count(p, 0). You can use degree of the polynomial for
 * an estimate of the number of roots.
 */
void lp_upolynomial_roots_isolate(const lp_upolynomial_t* p, lp_algebraic_number_t* roots, size_t* roots_size);

/**
 * Finds the root for an univariate polynomial in Zp. The array roots will be allocated,
 * and the user should de-allocate it. The size parameter will be updated with the size
 * of the array.
 */
void lp_upolynomial_roots_find_Zp(const lp_upolynomial_t* f, lp_integer_t** roots, size_t* roots_size);

/**
 * Reverses the coefficient of p in place. The result polynomial is
 * p'(x) = x^n * p(1/x) = a_n + ... + a_0 * x^n
 */
void lp_upolynomial_reverse_in_place(lp_upolynomial_t* p);


#ifdef __cplusplus
} /* close extern "C" { */
#endif
