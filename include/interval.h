/*
 * interval.h
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"
#include "value.h"

/**
 * An interval (a, b) with both point being values. This side is open is _open
 * is true.
 */
struct lp_interval_struct {
  /** Is the end at the point a open */
  size_t a_open : 1;
  /** Is the end at the point b open */
  size_t b_open : 1;
  /** Is this interval a point */
  size_t is_point : 1;
  /** The left end */
  lp_value_t a;
  /** The right end */
  lp_value_t b;
};

/** Construct the interval (a, b) */
void lp_interval_construct(lp_interval_t* I, const lp_value_t* a, int a_open, const lp_value_t* b, int b_open);

/** Construct the interval [0,0] */
void lp_interval_construct_zero(lp_interval_t* I);

/** Construct the interval [a, a] */
void lp_interval_construct_point(lp_interval_t* I, const lp_value_t* a);

/** Construct the interval (a, b) */
void lp_interval_construct_copy(lp_interval_t* I, const lp_interval_t* from);

/** Construct the interval (-inf, +inf) */
void lp_interval_construct_full(lp_interval_t* I);

/** Assign from another interval */
void lp_interval_assign(lp_interval_t* I, const lp_interval_t* from);

/** Set the lower bound */
void lp_interval_set_a(lp_interval_t* I, const lp_value_t* a, int a_open);

/** Set the upper bound */
void lp_interval_set_b(lp_interval_t* I, const lp_value_t* b, int b_open);

/** Collapse to a point v in the interval */
void lp_interval_collapse_to(lp_interval_t* I, const lp_value_t* v);

/** Destroy the interval */
void lp_interval_destruct(lp_interval_t* I);

/** Swap the two intervals */
void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2);

/** Check if the value is contained in the interval */
int lp_interval_contains(const lp_interval_t* I, const lp_value_t* v);

/** Returns an approximation of the log interval size */
int lp_interval_size_approx(const lp_interval_t* I);

/** Prints the interval to the given stream. */
int lp_interval_print(const lp_interval_t* I, FILE* out);

/** Returns the string representation of the interval */
char* lp_interval_to_string(const lp_interval_t* I);

/** Is this interval a point */
int lp_interval_is_point(const lp_interval_t* I);

/** Get the point value (it has to be a point) */
const lp_value_t* lp_interval_get_point(const lp_interval_t* I);

/** Get the upper bound */
const lp_value_t* lp_interval_get_lower_bound(const lp_interval_t* I);

/** Get the lower bound strict (open)? */
const lp_value_t* lp_interval_get_upper_bound(const lp_interval_t* I);

/** Returns a value in the interval */
void lp_interval_pick_value(const lp_interval_t* I, lp_value_t* v);

/** Compares the lower bounds of the intervals */
int lp_interval_cmp_lower_bounds(const lp_interval_t* I1, const lp_interval_t* I2);

/** Compares the uppoer bounds of the intervals */
int lp_interval_cmp_upper_bounds(const lp_interval_t* I1, const lp_interval_t* I2);

/**
 * Comparison of intervals based on upper bounds, with additional intersect info.
 */
typedef enum {
  /* I1: (  )
   * I2:      (   ) */
  LP_INTERVAL_CMP_LT_NO_INTERSECT,
  /* I1: (   )
   * I2:   (   )    */
  LP_INTERVAL_CMP_LT_WITH_INTERSECT,
  /* I1: (   )
   * I2: (     )    */
  LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1,
  /* I1: (     ]
   * I2:   (   ]    */
  LP_INTERVAL_CMP_LEQ_WITH_INTERSECT_I2,
  /* I1: (   ]
   * I2: (   ]      */
  LP_INTERVAL_CMP_EQ,
  /* I1:   (   ]
   * I2: (     ]    */
  LP_INTERVAL_CMP_GEQ_WITH_INTERSECT_I1,
  /* I1: (       )
   * I2: (    )     */
  LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2,
  /* I1:   (    )
   * I2: (    )     */
  LP_INTERVAL_CMP_GT_WITH_INTERSECT,
  /* I1:      (   )
   * I2: (  )       */
  LP_INTERVAL_CMP_GT_NO_INTERSECT
} lp_interval_cmp_t;

/**
 * Compares the two intervals.
 */
lp_interval_cmp_t lp_interval_cmp(const lp_interval_t* I1, const lp_interval_t* I2);

/**
 * Compares the two intervals and assigns the intersect (if any).
 */
lp_interval_cmp_t lp_interval_cmp_with_intersect(const lp_interval_t* I1, const lp_interval_t* I2, lp_interval_t* P);

