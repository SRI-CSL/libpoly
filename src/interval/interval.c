/*
 * dyadic_interval.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include "interval/interval.h"

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

const interval_ops_t interval_ops = {
    interval_construct,
    interval_construct_point,
    interval_construct_copy,
    interval_construct_from_dyadic,
    interval_construct_from_int,
    interval_construct_from_integer,
    interval_destruct,
    interval_swap,
    interval_print
};

const dyadic_interval_ops_t dyadic_interval_ops = {
    dyadic_interval_construct,
    dyadic_interval_construct_copy,
    dyadic_interval_construct_from_int,
    dyadic_interval_construct_from_integer,
    dyadic_interval_construct_point,
    dyadic_interval_construct_from_split,
    dyadic_interval_construct_intersection,
    dyadic_interval_destruct,
    dyadic_interval_swap,
    dyadic_interval_collapse_to,
    dyadic_interval_set_a,
    dyadic_interval_set_b,
    dyadic_interval_print,
    dyadic_interval_equals,
    dyadic_interval_contains,
    dyadic_interval_disjunct,
    dyadic_interval_scale
};
