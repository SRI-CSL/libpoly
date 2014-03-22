/*
 * interval.h
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#pragma once

#include "number/rational.h"
#include "number/dyadic_rational.h"

/**
 * An interval (a, b) with both points rationals. The side is open if
 * the _open is true.
 */
typedef struct {
  /** Is the end at the point a open */
  size_t a_open : 1;
  /** Is the end at the point b open */
  size_t b_open : 1;
  /** Is this a single point */
  size_t is_point : 1;
  /** The left end */
  rational_t a;
  /** The right end */
  rational_t b;
} interval_t;

typedef struct {

  /** Construct the interval (a, b) */
  void (*construct) (interval_t* I, const rational_t* a, int a_open, const rational_t* b, int b_open);

  /** Construct the interval [a, a] */
  void (*consturct_point) (interval_t* I, const rational_t* a);

  /** Construct the interval (a, b) */
  void (*construct_copy) (interval_t* I, const interval_t* from);

  /** Construct the interval (a, b) */
  void (*construct_from_dyadic) (interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_int) (interval_t* I, long a, int a_open, long b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_integer) (interval_t* I, const integer_t* a, int a_open, const integer_t* b, int b_open);

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
  dyadic_rational_t a;
  /** The right end */
  dyadic_rational_t b;
} dyadic_interval_t;

typedef struct {

  /** Construct the interval (a, b) */
  void (*construct) (dyadic_interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open);

  /** Construct the interval (a, b) */
  void (*construct_copy) (dyadic_interval_t* I, const dyadic_interval_t* from);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_int) (dyadic_interval_t* I, long a, int a_open, long b, int b_open);

  /** Construct from interval (a, b) from integers */
  void (*construct_from_integer) (dyadic_interval_t* I, const integer_t* a, int a_open, const integer_t* b, int b_open);

  /** Construct the interval [q, q] */
  void (*construct_point) (dyadic_interval_t* I, const dyadic_rational_t* q);

  /**
   * Construct the split of I into left and right in the middle. The mid point
   * will be open in left/right depending on the passed flags.
   */
  void (*construct_from_split) (dyadic_interval_t* I_left, dyadic_interval_t* I_right, const dyadic_interval_t* I, int left_open, int right_open);

  /** Construct an intersection of two non-disjoint intervals of the SAME SIZE. */
  void (*construct_intersection) (dyadic_interval_t* I, const dyadic_interval_t* I1, const dyadic_interval_t* I2);

  /** Destroy the interval */
  void (*destruct) (dyadic_interval_t* I);

  /** Swap two intervals */
  void (*swap) (dyadic_interval_t* I1, dyadic_interval_t* I2);

  /** Collapse the interval to a single point */
  void (*collapse_to) (dyadic_interval_t* I, const dyadic_rational_t* q);

  /** Set the left point of the interval */
  void (*set_a) (dyadic_interval_t* I, const dyadic_rational_t* a, int a_open);

  void (*set_b) (dyadic_interval_t* I, const dyadic_rational_t* b, int b_open);
  /** Set the right point of the interval */

  /** Prints the interval to the given stream. */
  int (*print) (const dyadic_interval_t* I, FILE* out);

  /** Check whether two intervals are equal  */
  int (*equals) (const dyadic_interval_t* I1, const dyadic_interval_t* I2);

  /** Check whether the number q is in the interval I */
  int (*contains) (const dyadic_interval_t* I, const dyadic_rational_t* q);

  /** Check whether two intervals are disjunct */
  int (*disjunct) (const dyadic_interval_t* I1, const dyadic_interval_t* I2);

  /** Scale the interval by power of 2 */
  void (*scale) (dyadic_interval_t* I, int n);

} dyadic_interval_ops_t;

extern const dyadic_interval_ops_t dyadic_interval_ops;
