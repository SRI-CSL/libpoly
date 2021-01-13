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

#include "dyadic_interval.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Algebraic number represented as the only root of the polynomial f in the
 * interval (a,b). The signs at the points a and b are kept to improve
 * refinement of the interval when needed. If f is 0, then the interval is
 * a single point, and that is the value of the number.
 */
struct lp_algebraic_number_struct {
  lp_upolynomial_t* f;
  lp_dyadic_interval_t I;
  int sgn_at_a, sgn_at_b;
};

/**
 * Construct the algebraic number given its polynomial and the isolating
 * interval. The number takes over the reference of f.
 */
void lp_algebraic_number_construct(lp_algebraic_number_t* a, lp_upolynomial_t* f, const lp_dyadic_interval_t* I);

/** Construct a zero algebraic number */
void lp_algebraic_number_construct_zero(lp_algebraic_number_t* a);

/** Construct a zero algebraic number */
void lp_algebraic_number_construct_one(lp_algebraic_number_t* a);

/** Construct a copy of the algebraic number. */
void lp_algebraic_number_construct_copy(lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

/** Construct the algebraic number from an integerl */
void lp_algebraic_number_construct_from_integer(lp_algebraic_number_t* a, const lp_integer_t* z);

/** Construct the algebraic number from a dyadic rational */
void lp_algebraic_number_construct_from_dyadic_rational(lp_algebraic_number_t* a, const lp_dyadic_rational_t* q);

/** Construct the algebraic number from a rational */
void lp_algebraic_number_construct_from_rational(lp_algebraic_number_t* a, const lp_rational_t* q);

/** Destruct the number */
void lp_algebraic_number_destruct(lp_algebraic_number_t* a);

/** Swap the two numbers */
void lp_algebraic_number_swap(lp_algebraic_number_t* a, lp_algebraic_number_t* b);

/** Get the sign of the algebraic number */
int lp_algebraic_number_sgn(const lp_algebraic_number_t* a);

/** Compare two algebraic numbers */
int lp_algebraic_number_cmp(const lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2);

/** Compare algebraic number to an integer */
int lp_algebraic_number_cmp_integer(const lp_algebraic_number_t* a1, const lp_integer_t* a2);

/** Compare algebraic number to a dyadic rational */
int lp_algebraic_number_cmp_dyadic_rational(const lp_algebraic_number_t* a1, const lp_dyadic_rational_t* a2);

/** Compare algebraic number to a rational */
int lp_algebraic_number_cmp_rational(const lp_algebraic_number_t* a1, const lp_rational_t* a2);

/** Void version of the comparison, use with care. */
int lp_algebraic_number_cmp_void(const void* a1, const void* a2);

/** Print the number */
int lp_algebraic_number_print(const lp_algebraic_number_t* a, FILE* out);

/** Return a string representation of the number */
char* lp_algebraic_number_to_string(const lp_algebraic_number_t* a);

/** Convert to double approximation */
double lp_algebraic_number_to_double(const lp_algebraic_number_t* a);

/** Convert to rational approximation */
void lp_algebraic_number_to_rational(const lp_algebraic_number_t* a, lp_rational_t* q);

/** Get the midpoint of the defining interval */
void lp_algebraic_number_get_dyadic_midpoint(const lp_algebraic_number_t* a, lp_dyadic_rational_t* q);

/** Get the midpoint of the defining interval */
void lp_algebraic_number_get_rational_midpoint(const lp_algebraic_number_t* a, lp_rational_t* q);

/** Refine the number by halving its interval. */
void lp_algebraic_number_refine(lp_algebraic_number_t* a);

/**
 * Same as above, but const version for convenience: NOT CONST, the number is
 * the same but internally it might change.
 */
void lp_algebraic_number_refine_const(const lp_algebraic_number_t* a);

/** Restore the interval that has been lost due to refinement (be careful) */
void lp_algebraic_number_restore_interval(lp_algebraic_number_t* a, const lp_dyadic_interval_t* I);

/** Same as above, but const */
void lp_algebraic_number_restore_interval_const(const lp_algebraic_number_t* a, const lp_dyadic_interval_t* I);


/** Addition */
void lp_algebraic_number_add(lp_algebraic_number_t* sum, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Subtraction */
void lp_algebraic_number_sub(lp_algebraic_number_t* sub, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Negation */
void lp_algebraic_number_neg(lp_algebraic_number_t* neg, const lp_algebraic_number_t* a);

/** Multiplication */
void lp_algebraic_number_mul(lp_algebraic_number_t* mul, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Multiplicative inverse 1/a (a != 0) */
void lp_algebraic_number_inv(lp_algebraic_number_t* inv, const lp_algebraic_number_t* a);

/** Division (b != 0) */
void lp_algebraic_number_div(lp_algebraic_number_t* div, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b);

/** Exponentiation */
void lp_algebraic_number_pow(lp_algebraic_number_t* pow, const lp_algebraic_number_t* a, unsigned n);

/** Returns true if a rational number (not complete) */
int lp_algebraic_number_is_rational(const lp_algebraic_number_t* a);

/** Returns true if an integer number (complete) */
int lp_algebraic_number_is_integer(const lp_algebraic_number_t* a);

/** Returns the ceiling of the number */
void lp_algebraic_number_ceiling(const lp_algebraic_number_t* a, lp_integer_t* a_ceiling);

/** Returns the floor of the number */
void lp_algebraic_number_floor(const lp_algebraic_number_t* a, lp_integer_t* a_floor);

/** Returns a hash of the of the dyadic approximation of a */
size_t lp_algebraic_number_hash_approx(const lp_algebraic_number_t* a, unsigned precision);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
