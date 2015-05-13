/*
 * interval.h
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "rational.h"
#include "dyadic_rational.h"

/**
 * An interval (a, b) with both points rationals. The side is open if
 * the _open is true.
 */
typedef struct lp_interval_struct {
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
} lp_interval_t;


/**
 * A dyadic interval is an interval either an interval of the form (a, b) where
 * both a and b are dyadic rationals. This is different from the traditional
 * use of "dyadic interval", where intervals are always of the form (a/2^n, (a+1)/2^n).
 * This is so that we can perform interval arithmetic with the intervals
 * (otherwise, multiplication with scalars would not be allowed.
 */
typedef struct {
  /** Is the end at the point a open */
  size_t a_open : 1;
  /** Is the end at the point b open */
  size_t b_open : 1;
  /** Is this interval a point */
  size_t is_point : 1;
  /** The left end */
  lp_dyadic_rational_t a;
  /** The right end */
  lp_dyadic_rational_t b;
} lp_dyadic_interval_t;


/** Construct the interval (a, b) */
void lp_interval_construct(lp_interval_t* I, const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open);

/** Construct the interval [0,0] */
void lp_interval_construct_zero(lp_interval_t* I);

/** Construct the interval [a, a] */
void lp_interval_construct_point(lp_interval_t* I, const lp_rational_t* a);

/** Construct the interval (a, b) */
void lp_interval_construct_copy(lp_interval_t* I, const lp_interval_t* from);

/** Construct the interval (a, b) */
void lp_interval_construct_from_dyadic(lp_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

/** Construct from diadic interval */
void lp_interval_construct_from_dyadic_interval(lp_interval_t* I, const lp_dyadic_interval_t* from);

/** Construct from interval (a, b) from integers */
void lp_interval_construct_from_int(lp_interval_t* I, long a, int a_open, long b, int b_open);

/** Construct from interval (a, b) from integers */
void lp_interval_construct_from_integer(lp_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

/** Assign from another interval */
void lp_interval_assign(lp_interval_t* I, const lp_interval_t* from);

/** Destroy the interval */
void lp_interval_destruct(lp_interval_t* I);

/** Swap the two intervals */
void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2);

/** Prints the interval to the given stream. */
int lp_interval_print(const lp_interval_t* I, FILE* out);

/** Check if the interval contains q */
int lp_interval_contains(const lp_interval_t* I, const lp_rational_t* q);

/** Check whether the interval contains 0 */
int lp_interval_contains_zero(const lp_interval_t* I);

/** Returns the "sign" of the interval: 0 in 0 in I, negative if I < 0, positive if I > 0. */
int lp_interval_sgn(const lp_interval_t* I);



/** Construct the interval (a, b) */
void lp_dyadic_interval_construct(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

/** Construct the interval [0, 0] */
void lp_dyadic_interval_construct_zero(lp_dyadic_interval_t* I);

/** Construct the interval (a, b) */
void lp_dyadic_interval_construct_copy(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from);

/** Construct from interval (a, b) from integers */
void lp_dyadic_interval_construct_from_int(lp_dyadic_interval_t* I, long a, int a_open, long b, int b_open);

/** Construct from interval (a, b) from integers */
void lp_dyadic_interval_construct_from_integer(lp_dyadic_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

/** Construct the interval [q, q] */
void lp_dyadic_interval_construct_point(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

/**
 * Construct the split of I into left and right in the middle. The mid point
 * will be open in left/right depending on the passed flags.
 */
void lp_dyadic_interval_construct_from_split(lp_dyadic_interval_t* I_left, lp_dyadic_interval_t* I_right, const lp_dyadic_interval_t* I, int left_open, int right_open);

/** Construct an intersection of two non-disjoint intervals of the SAME SIZE. */
void lp_dyadic_interval_construct_intersection(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

/** Assign from another interval */
void lp_dyadic_interval_assign(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from);

/** Destroy the interval */
void lp_dyadic_interval_destruct(lp_dyadic_interval_t* I);

/** Swap two intervals */
void lp_dyadic_interval_swap(lp_dyadic_interval_t* I1, lp_dyadic_interval_t* I2);

/** Collapse the interval to a single point */
void lp_dyadic_interval_collapse_to(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

/** Set the left point of the interval */
void lp_dyadic_interval_set_a(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open);

/** Set the right point of the interval */
void lp_dyadic_interval_set_b(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* b, int b_open);

/** Prints the interval to the given stream. */
int lp_dyadic_interval_print(const lp_dyadic_interval_t* I, FILE* out);

/** Check whether two intervals are equal  */
int lp_dyadic_interval_equals(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

/** Check whether the number q is in the interval I */
int lp_dyadic_interval_contains(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

/** Check whether the interval contains 0 */
int lp_dyadic_interval_contains_zero(const lp_dyadic_interval_t* I);

/** Returns the "sign" of the interval: 0 in 0 in I, negative if I < 0, positive if I > 0. */
int lp_dyadic_interval_sgn(const lp_dyadic_interval_t* I);

/**
 * Compare the interval to a integer, returns 0 if number is inside,
 * -1 if interval is below number, or +1 if interval is above number.
 */
int lp_dyadic_interval_cmp_integer(const lp_dyadic_interval_t* I, const lp_integer_t* z);

/**
 * Compare the interval to a dyadic rational, returns 0 if number is inside,
 * -1 if interval is below number, or +1 if interval is above number.
 */
int lp_dyadic_interval_cmp_dyadic_rational(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

/**
 * Compare the interval to a rational, returns 0 if number is inside,
 * -1 if interval is below number, or +1 if interval is above number.
 */
int lp_dyadic_interval_cmp_rational(const lp_dyadic_interval_t* I, const lp_rational_t* q);


/** Check whether two intervals are disjunct */
int lp_dyadic_interval_disjunct(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

/** Scale the interval by power of 2 */
void lp_dyadic_interval_scale(lp_dyadic_interval_t* I, int n);

/** Is this interval a point */
int lp_dyadic_interval_is_point(const lp_dyadic_interval_t* I);

/** Get the point value */
const lp_dyadic_rational_t* lp_dyadic_interval_get_point(const lp_dyadic_interval_t* I);

/** Is this interval a point */
int lp_interval_is_point(const lp_interval_t* I);

/** Get the point value */
const lp_rational_t* lp_interval_get_point(const lp_interval_t* I);

