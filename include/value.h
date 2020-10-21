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
#include "dyadic_rational.h"
#include "algebraic_number.h"
#include "value.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Types of values for the assignment */
typedef enum {
  /** No value */
  LP_VALUE_NONE,
  /** Integer value */
  LP_VALUE_INTEGER,
  /** Dyadic rational */
  LP_VALUE_DYADIC_RATIONAL,
  /** Rational number */
  LP_VALUE_RATIONAL,
  /** Reduced algebraic number (univariate representation) */
  LP_VALUE_ALGEBRAIC,
  /** Special value +inf */
  LP_VALUE_PLUS_INFINITY,
  /** Special value -inf */
  LP_VALUE_MINUS_INFINITY,
} lp_value_type_t;

/** A value is a choice of the available types */
typedef union {
  lp_integer_t z;
  lp_rational_t q;
  lp_dyadic_rational_t dy_q;
  lp_algebraic_number_t a;
} lp_value_union_t;

/** A value is a tagged union of available type */
struct lp_value_struct {
  lp_value_type_t type;
  lp_value_union_t value;
};

/** Construct a value */
void lp_value_construct(lp_value_t* v, lp_value_type_t type, const void* data);

/** Construct a zero value */
void lp_value_construct_zero(lp_value_t* v);

/** Consruct and integer */
void lp_value_construct_int(lp_value_t* v, long x);

/** Construct the null value */
void lp_value_construct_none(lp_value_t* v);

/** Returns a null value */
const lp_value_t* lp_value_none(void);

/** Returns -inf */
const lp_value_t* lp_value_minus_infinity(void);

/** Returns +inf */
const lp_value_t* lp_value_plus_infinity(void);

/** Construct a copy of the given value */
void lp_value_construct_copy(lp_value_t* v, const lp_value_t* from);

/** Allocate and construct */
lp_value_t* lp_value_new(lp_value_type_t type, const void* data);

/** Allocate and construc a copy */
lp_value_t* lp_value_new_copy(const lp_value_t* from);

/** Destruct the value */
void lp_value_destruct(lp_value_t* v);

/** Destruct and free the value */
void lp_value_delete(lp_value_t* v);

/** Get a hash of the value (not a good hash == lp_value_hash_approx(v, 0) */
size_t lp_value_hash(const lp_value_t* v);

/** Returns a hash of the of the dyadic approximation of v */
size_t lp_value_hash_approx(const lp_value_t* v, unsigned precision);

/** Assign */
void lp_value_assign(lp_value_t* v, const lp_value_t* from);

/** Assign 0 */
void lp_value_assign_zero(lp_value_t* v);

/** Assign a to value */
void lp_value_assign_raw(lp_value_t* v, lp_value_type_t type, const void* data);

/** Swap two values */
void lp_value_swap(lp_value_t* v1, lp_value_t* v2);

/** Get the approximate value */
void lp_value_approximate(const lp_value_t* v, lp_rational_interval_t* approx);

/** Compare two values. */
int lp_value_cmp(const lp_value_t* v1, const lp_value_t* v2);

/** Void version of the comparison, use with care. */
int lp_value_cmp_void(const void* v1, const void* v2);

/** Compare to a rational */
int lp_value_cmp_rational(const lp_value_t* v, const lp_rational_t* q);

/** Print the value */
int lp_value_print(const lp_value_t* v, FILE* out);

/** Return a string representation */
char* lp_value_to_string(const lp_value_t* v);

/** Sign of the value */
int lp_value_sgn(const lp_value_t* v);


/**
 * Check if the value is a rational number: either an integer, dyadic rational,
 * a rational, or a algebraic number that has reduced to a point.
 */
int lp_value_is_rational(const lp_value_t* v);

/**
 * Check if the value is an integer number: either an integer, dyadic rational,
 * a rational, or a algebraic number that has reduced to a point.
 */
int lp_value_is_integer(const lp_value_t* v);


/**
 * Check if the value is +/- infinity.
 */
int lp_value_is_infinity(const lp_value_t* v);

/** Get the floor of the value */
void lp_value_floor(const lp_value_t* v, lp_integer_t* v_floor);

/** Get the ceiling of the value */
void lp_value_ceiling(const lp_value_t* v, lp_integer_t* v_ceiling);

/**
 * Get the rational if is_rational is true.
 */
void lp_value_get_rational(const lp_value_t* v, lp_rational_t* q);

/** Get the numerator (only if integer, dyadic rational, or rational) */
void lp_value_get_num(const lp_value_t* v, lp_integer_t* num);

/** Get the denominator (only if integer, dyadic rational, or rational) */
void lp_value_get_den(const lp_value_t* v, lp_integer_t* den);

/** Get a value in the interval [a, b] with strictness of the interval given by a_strict and b_strict. */
void lp_value_get_value_between(const lp_value_t* a, int a_strict, const lp_value_t* b, int b_strict, lp_value_t* v);

/** Get an approximation of the size between lower and upper */
int lp_value_get_distance_size_approx(const lp_value_t* lower, const lp_value_t* upper);

/** Get the double (approximation) of the value */
double lp_value_to_double(const lp_value_t* v);

/** Addition. Does not support adding -inf and +inf. */
void lp_value_add(lp_value_t* sum, const lp_value_t* a, const lp_value_t* b);

/** Subtraction. Does not support subtracting inf and inf of the same sign */
void lp_value_sub(lp_value_t* sub, const lp_value_t* a, const lp_value_t* b);

/** Negation */
void lp_value_neg(lp_value_t* neg, const lp_value_t* a);

/** Multiplication. Does not support multiplying 0 and inf. */
void lp_value_mul(lp_value_t* mul, const lp_value_t* a, const lp_value_t* b);

/** Multiplicative inverse. Does not support ./0, inf/inf. */
void lp_value_inv(lp_value_t* inv, const lp_value_t* a);

/** Division (b != 0) */
void lp_value_div(lp_value_t* div, const lp_value_t* a, const lp_value_t* b);

/** Exponentiation */
void lp_value_pow(lp_value_t* pow, const lp_value_t* a, unsigned n);


#ifdef __cplusplus
} /* close extern "C" { */
#endif
