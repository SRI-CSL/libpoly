/*
 * value.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "integer.h"
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"
#include "value.h"

#include <stdio.h>

/** Types of values for the assignment */
typedef enum {
  /** No value */
  LP_VALUE_NONE,
  /** Integer value */
  LP_VALUE_INTEGER,
  /** Dyadic rational */
  LP_VALUE_DYADIC_RATIONAL,
  /** Ratioanl number */
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

/** Construct the null value */
void lp_value_construct_none(lp_value_t* v);

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

/** Assign */
void lp_value_assign(lp_value_t* v, const lp_value_t* from);

/** Assign a to value */
void lp_value_assign_raw(lp_value_t* v, lp_value_type_t type, const void* data);

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

/**
 * Check if the value is a rational number: either an integer, dyadic rational,
 * a rational, or a algebraic number that has reduced to a point.
 */
int lp_value_is_rational(const lp_value_t* v);

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

