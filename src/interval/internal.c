/*
 * interval_internal.c
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#include "interval/internal.h"

#include <assert.h>
#include <limits.h>

/** Prints the interval to the given stream. */
int interval_print(const interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += rational_ops.print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += rational_ops.print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += rational_ops.print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

int dyadic_interval_print(const dyadic_interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += dyadic_rational_ops.print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += dyadic_rational_ops.print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += dyadic_rational_ops.print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}
