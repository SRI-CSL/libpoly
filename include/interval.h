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
typedef struct interval_struct {
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
} interval_t;

typedef struct {

  /** Construct the interval (a, b) */
  void (*construct) (interval_t* I, const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open);

  /** Construct the interval [a, a] */
  void (*consturct_point) (interval_t* I, const lp_rational_t* a);

  /** Construct the interval (a, b) */
  void (*construct_copy) (interval_t* I, const interval_t* from);

  /** Construct the interval (a, b) */
  void (*construct_from_dyadic) (interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_int) (interval_t* I, long a, int a_open, long b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_integer) (interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

  /** Destroy the interval */
  void (*destruct) (interval_t* I);

  /** Swap the two intervals */
  void (*swap) (interval_t* I1, interval_t* I2);

  /** Prints the interval to the given stream. */
  int (*print) (const interval_t* I, FILE* out);

} interval_ops_t;

extern const interval_ops_t interval_ops;

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

typedef struct {

  /** Construct the interval (a, b) */
  void (*construct) (lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open);

  /** Construct the interval (a, b) */
  void (*construct_copy) (lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_int) (lp_dyadic_interval_t* I, long a, int a_open, long b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_integer) (lp_dyadic_interval_t* I, const lp_integer_t* a, int a_open, const lp_integer_t* b, int b_open);

  /** Construct the interval [q, q] */
  void (*construct_point) (lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

  /**
   * Construct the split of I into left and right in the middle. The mid point
   * will be open in left/right depending on the passed flags.
   */
  void (*construct_from_split) (lp_dyadic_interval_t* I_left, lp_dyadic_interval_t* I_right, const lp_dyadic_interval_t* I, int left_open, int right_open);

  /** Construct an intersection of two non-disjoint intervals of the SAME SIZE. */
  void (*construct_intersection) (lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

  /** Destroy the interval */
  void (*destruct) (lp_dyadic_interval_t* I);

  /** Swap two intervals */
  void (*swap) (lp_dyadic_interval_t* I1, lp_dyadic_interval_t* I2);

  /** Collapse the interval to a single point */
  void (*collapse_to) (lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

  /** Set the left point of the interval */
  void (*set_a) (lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open);

  void (*set_b) (lp_dyadic_interval_t* I, const lp_dyadic_rational_t* b, int b_open);
  /** Set the right point of the interval */

  /** Prints the interval to the given stream. */
  int (*print) (const lp_dyadic_interval_t* I, FILE* out);

  /** Check whether two intervals are equal  */
  int (*equals) (const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

  /** Check whether the number q is in the interval I */
  int (*contains) (const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q);

  /** Check whether two intervals are disjunct */
  int (*disjunct) (const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2);

  /** Scale the interval by power of 2 */
  void (*scale) (lp_dyadic_interval_t* I, int n);

} lp_dyadic_interval_ops_t;

extern const lp_dyadic_interval_ops_t lp_dyadic_interval_ops;
