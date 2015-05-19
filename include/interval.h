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

/** Destroy the interval */
void lp_interval_destruct(lp_interval_t* I);

/** Swap the two intervals */
void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2);

/** Prints the interval to the given stream. */
int lp_interval_print(const lp_interval_t* I, FILE* out);

/** Returns the string representation of the interval */
char* lp_interval_to_string(const lp_interval_t* I);

/** Is this interval a point */
int lp_interval_is_point(const lp_interval_t* I);

/** Get the point value */
const lp_value_t* lp_interval_get_point(const lp_interval_t* I);
