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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Sign condition contains two signs, both either -1, 0, 1, in ascending order.
 * If the second sign is equal to the first one, then only the first one is
 * relevant.
 */
enum lp_sign_condition_enum {
  /** <  0 */
  LP_SGN_LT_0,
  /** <= 0 */
  LP_SGN_LE_0,
  /** == 0 */
  LP_SGN_EQ_0,
  /** != 0 */
  LP_SGN_NE_0,
  /** >  0 */
  LP_SGN_GT_0,
  /** >= 0 */
  LP_SGN_GE_0
};

typedef enum lp_sign_condition_enum lp_sign_condition_t;

/**
 * Check if the sign condition is consistent with the given sign.
 */
int lp_sign_condition_consistent(lp_sign_condition_t sgn_condition, int sign);

/**
 * CHeck if the sign condition is consistent with the given interval.
 */
int lp_sign_condition_consistent_interval(lp_sign_condition_t sgn_condition, const lp_interval_t* I);

/**
 * Negate the sign condition.
 */
lp_sign_condition_t lp_sign_condition_negate(lp_sign_condition_t sgn_condition);

/**
 * Print the sign condition.
 */
int lp_sign_condition_print(lp_sign_condition_t sgn_condition, FILE* out);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
