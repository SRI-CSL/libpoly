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

#include <sign_condition.h>

lp_sign_condition_t lp_sign_condition_negate(lp_sign_condition_t sgn_condition) {
  switch (sgn_condition) {
  case LP_SGN_LT_0:
    return LP_SGN_GE_0;
  case LP_SGN_LE_0:
    return LP_SGN_GT_0;
  case LP_SGN_EQ_0:
    return LP_SGN_NE_0;
  case LP_SGN_NE_0:
    return LP_SGN_EQ_0;
  case LP_SGN_GT_0:
    return LP_SGN_LE_0;
  case LP_SGN_GE_0:
    return LP_SGN_LT_0;
  }
  // Just for compilers
  return LP_SGN_EQ_0;
}

int lp_sign_condition_print(lp_sign_condition_t sgn_condition, FILE* out) {
  int ret = 0;
  switch (sgn_condition) {
  case LP_SGN_LT_0:
    ret += fprintf(out, "< 0");
    break;
  case LP_SGN_LE_0:
    ret += fprintf(out, "<= 0");
    break;
  case LP_SGN_EQ_0:
    ret += fprintf(out, "== 0");
    break;
  case LP_SGN_NE_0:
    ret += fprintf(out, "!= 0");
    break;
  case LP_SGN_GT_0:
    ret += fprintf(out, "> 0");
    break;
  case LP_SGN_GE_0:
    ret += fprintf(out, ">= 0");
    break;
  }
  return ret;
}

int lp_sign_condition_consistent(lp_sign_condition_t sgn_condition, int sign) {
  switch (sgn_condition) {
  case LP_SGN_LT_0:
    return sign < 0;
  case LP_SGN_LE_0:
    return sign <= 0;
  case LP_SGN_EQ_0:
    return sign == 0;
  case LP_SGN_NE_0:
    return sign != 0;
  case LP_SGN_GT_0:
    return sign > 0;
  case LP_SGN_GE_0:
    return sign >= 0;
  }
  return 0;
}

int lp_sign_condition_consistent_interval(lp_sign_condition_t sgn_condition, const lp_interval_t* I) {
  (void) I;
  (void) sgn_condition;
  return 0;
}

