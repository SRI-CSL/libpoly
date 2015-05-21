/*
 * feasibility_set.c
 *
 *  Created on: Feb 24, 2015
 *      Author: dejan
 */

#include <feasibility_set.h>
#include <stdlib.h>
#include <assert.h>

#include "poly.h"
#include "value.h"
#include "interval.h"

#include "polynomial/polynomial.h"
#include "utils/debug_trace.h"
#include "polynomial/feasibility_set.h"

void lp_feasibility_set_construct(lp_feasibility_set_t* s, size_t size) {
  s->size = 0;
  s->capacity = size;
  if (size) {
    s->intervals = malloc(s->capacity * sizeof(lp_interval_t));
  } else {
    s->intervals = 0;
  }
}

lp_feasibility_set_t* lp_feasibility_set_new_internal(size_t size) {
  lp_feasibility_set_t* result = malloc(sizeof(lp_feasibility_set_t));
  lp_feasibility_set_construct(result, size);
  return result;
}

void lp_feasibility_set_destruct(lp_feasibility_set_t* s) {
  size_t i;
  for (i = 0; i < s->size; ++ i) {
    lp_interval_destruct(s->intervals + i);
  }
  free(s->intervals);
}

lp_feasibility_set_t* lp_feasibility_set_new() {
  lp_feasibility_set_t* result = lp_feasibility_set_new_internal(1);
  lp_value_t inf_neg, inf_pos;
  lp_value_construct(&inf_neg, LP_VALUE_MINUS_INFINITY, 0);
  lp_value_construct(&inf_pos, LP_VALUE_PLUS_INFINITY, 0);
  lp_interval_construct(result->intervals, &inf_neg, 1, &inf_pos, 1);
  result->size = 1;
  return result;
}

void lp_feasibility_set_delete(lp_feasibility_set_t* set) {
  lp_feasibility_set_destruct(set);
  free(set);
}

int lp_feasibility_set_is_empty(const lp_feasibility_set_t* set) {
  return set->size == 0;
}

int lp_feasibility_set_print(const lp_feasibility_set_t* set, FILE* out) {
  int ret = 0;
  size_t i;
  ret += fprintf(out, "{ ");
  for(i = 0; i < set->size; ++ i) {
    if (i) {
      ret += fprintf(out, ", ");
    }
    ret += lp_interval_print(&set->intervals[i], out);
  }
  ret += fprintf(out, " }");
  return ret;
}

int lp_feasibility_set_contains(const lp_feasibility_set_t* set, const lp_value_t* value) {
  // TODO: binary search
  for (size_t i = 0; i < set->size; ++ i) {
    if (lp_interval_contains(set->intervals + i, value)) {
      return 1;
    }
  }
  return 0;
}

lp_value_t* lp_feasibility_set_pick_value(const lp_feasibility_set_t* set) {
  (void)set;
  assert(0);
  return 0;
}
