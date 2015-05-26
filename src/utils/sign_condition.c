/*
 * sign_condition.c
 *
 *  Created on: May 26, 2015
 *      Author: dejan
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
    ret += fprintf(out, "= 0");
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
