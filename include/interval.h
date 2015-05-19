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

