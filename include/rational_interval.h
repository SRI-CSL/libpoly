/*
 * rational_interval.h
 *
 *  Created on: May 19, 2015
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "rational.h"

/**
 * An interval (a, b) with both points rationals. The side is open if
 * the _open is true.
 */
typedef struct lp_rational_interval_struct {
  /** Is the end at the point a open */
  size_t a_open : 1;
  /** Is the end at the point b open */
  size_t b_open : 1;
  /** Is this a single point */
  size_t is_point : 1;
  /** The left end */
  lp_rational_t a;
  /** The right end */
  lp_rational_t b;
} lp_rational_interval_t;

/** Construct the interval (a, b) */
void lp_rational_interval_construct(lp_rational_interval_t* I, const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open);

/** Construct the interval [0,0] */
void lp_rational_interval_construct_zero(lp_rational_interval_t* I);

/** Construct the interval [a, a] */
void lp_rational_interval_construct_point(lp_rational_interval_t* I, const lp_rational_t* a);

/** Construct the interval (a, b) */
void lp_rational_interval_construct_copy(lp_rational_interval_t* I, const lp_rational_interval_t* from);

/** Construct the interval (a, b) */
void lp_rational_interval_construct_from_dyadic(lp_rational_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

/** Construct from diadic interval */
void lp_rational_interval_construct_from_dyadic_interval(lp_rational_interval_t* I, const lp_dyadic_interval_t* from);

/** Construct from interval (a, b) from integers */
void lp_rational_interval_construct_from_int(lp_rational_interval_t* I, long a, int a_open, long b, int b_open);

/** Construct from interval (a, b) from integers */
void lp_rational_interval_construct_from_integer(lp_rational_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

/** Assign from another interval */
void lp_rational_interval_assign(lp_rational_interval_t* I, const lp_rational_interval_t* from);

/** Destroy the interval */
void lp_rational_interval_destruct(lp_rational_interval_t* I);

/** Swap the two intervals */
void lp_rational_interval_swap(lp_rational_interval_t* I1, lp_rational_interval_t* I2);

/** Prints the interval to the given stream. */
int lp_rational_interval_print(const lp_rational_interval_t* I, FILE* out);

/** Check if the interval contains integer */
int lp_rational_interval_contains_integer(const lp_rational_interval_t* I, const lp_integer_t* z);

/** Check if the interval contains dyadic rational q */
int lp_rational_interval_contains_dyadic_rational(const lp_rational_interval_t* I, const lp_dyadic_rational_t* dy_q);

/** Check if the interval contains rational q */
int lp_rational_interval_contains_rational(const lp_rational_interval_t* I, const lp_rational_t* q);

/** Check if the interval contains the given algebraic */
int lp_rational_interval_contains_algebraic_number(const lp_rational_interval_t* I, const lp_algebraic_number_t* a);

/** Check if the interval contains the given value */
int lp_rational_interval_contains_value(const lp_rational_interval_t* I, const lp_value_t* v);

/** Check whether the interval contains 0 */
int lp_rational_interval_contains_zero(const lp_rational_interval_t* I);

/** Returns the "sign" of the interval: 0 in 0 in I, negative if I < 0, positive if I > 0. */
int lp_rational_interval_sgn(const lp_rational_interval_t* I);

/** Is this interval a point */
int lp_rational_interval_is_point(const lp_rational_interval_t* I);

/** Get the point value */
const lp_rational_t* lp_rational_interval_get_point(const lp_rational_interval_t* I);
