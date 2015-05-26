/*
 * feasibility_set.c
 *
 *  Created on: Feb 24, 2015
 *      Author: dejan
 */

#include <feasibility_set.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

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

lp_feasibility_set_t* lp_feasibility_set_new_from_intervals(lp_interval_t* intervals, size_t intervals_size) {
  lp_feasibility_set_t* result = lp_feasibility_set_new_internal(intervals_size);
  memcpy(result->intervals, intervals, sizeof(lp_interval_t)*intervals_size);
  result->size = intervals_size;
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

void lp_feasibility_set_construct_copy(lp_feasibility_set_t* set, const lp_feasibility_set_t* from) {
  lp_feasibility_set_construct(set, from->size);
  size_t i;
  for (i = 0; i < from->size; ++ i) {
    lp_interval_construct_copy(set->intervals + i, from->intervals + i);
  }
  set->size = from->size;
}

lp_feasibility_set_t* lp_feasibility_set_new_copy(const lp_feasibility_set_t* set) {
  lp_feasibility_set_t* result = malloc(sizeof(lp_feasibility_set_t));
  lp_feasibility_set_construct_copy(result, set);
  return result;
}

void lp_feasibiliy_set_assign(lp_feasibility_set_t* set, const lp_feasibility_set_t* from) {
  if (set != from) {
    lp_feasibility_set_destruct(set);
    lp_feasibility_set_construct_copy(set, from);
  }
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

char* lp_feasibility_set_to_string(const lp_feasibility_set_t* set) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_feasibility_set_print(set, f);
  fclose(f);
  return str;
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

void lp_feasibility_set_pick_value(const lp_feasibility_set_t* set, lp_value_t* value) {
  (void)set;
  (void)value;
  assert(0);
}

lp_feasibility_set_t* lp_feasibility_set_intersect(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2) {
  lp_feasibility_set_intersect_status_t status;
  return lp_feasibility_set_intersect_with_status(s1, s2, &status);
}

lp_feasibility_set_t* lp_feasibility_set_intersect_with_status(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2, lp_feasibility_set_intersect_status_t* status) {

  // Size of the result is at most max of the sizes
  size_t intervals_capacity = s1->size > s2->size ? s1->size : s2->size;
  // one extra for the working copy
  intervals_capacity ++;
  lp_interval_t* intervals = malloc(sizeof(lp_interval_t)*intervals_capacity);
  size_t intervals_size = 0;

  // Construct the intervals
  size_t i;
  for (i = 0; i < intervals_capacity; ++ i) {
    lp_interval_construct_zero(intervals + i);
  }

  int all_s1 = 1; // Result is s1
  int all_s2 = 1; // Result is s2

  // Scan from left to right and construct
  size_t s1_i = 0, s2_i = 0;
  for (; s1_i < s1->size && s2_i < s2->size;) {

    assert(intervals_size < intervals_capacity);

    lp_interval_t* P = intervals + intervals_size;

    if (trace_is_enabled("feasibility_set")) {
      tracef("s1[%zu] = ", s1_i); lp_interval_print(s1->intervals + s1_i, trace_out); tracef("\n");
      tracef("s2[%zu] = ", s2_i); lp_interval_print(s2->intervals + s2_i, trace_out); tracef("\n");
    }

    // Compare the current intervals
    lp_interval_cmp_t cmp = lp_interval_cmp_with_intersect(s1->intervals + s1_i, s2->intervals + s2_i, P);

    if (trace_is_enabled("feasibility_set")) {
      switch (cmp) {
      case LP_INTERVAL_CMP_LT_NO_INTERSECT:
      case LP_INTERVAL_CMP_GT_NO_INTERSECT:
        tracef("no intersect\n");
        break;
      default:
        tracef("intersect P = "); lp_interval_print(P, trace_out); tracef("\n");
        break;
      }
    }

    // Advance in at least one interval and take in the intersect if any
    switch (cmp) {
    case LP_INTERVAL_CMP_LT_NO_INTERSECT:
      /* I1: (  )
       * I2:      (   ) */
      s1_i ++;
      all_s1 = 0;
      break;
    case LP_INTERVAL_CMP_LT_WITH_INTERSECT:
      /* I1: (   )
       * I2:   (   )    */
      s1_i ++;
      all_s1 = 0;
      all_s2 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1:
      /* I1: (   )
       * I2: (     )    */
      s1_i ++;
      all_s2 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_LEQ_WITH_INTERSECT_I2:
      /* I1: (     ]
       * I2:   (   ]    */
      s1_i ++;
      s2_i ++;
      all_s1 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_EQ:
      /* I1: (   ]
       * I2: (   ]      */
      s1_i ++;
      s2_i ++;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_GEQ_WITH_INTERSECT_I1:
      /* I1:   (   ]
       * I2: (     ]    */
      s1_i ++;
      s2_i ++;
      all_s2 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2:
      /* I1: (       )
       * I2: (    )     */
      s2_i ++;
      all_s1 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_GT_WITH_INTERSECT:
      /* I1:   (    )
       * I2: (    )     */
      s2_i ++;
      all_s1 = 0;
      all_s2 = 0;
      intervals_size ++;
      break;
    case LP_INTERVAL_CMP_GT_NO_INTERSECT:
      /* I1:      (   )
       * I2: (  )       */
      s2_i ++;
      all_s2 = 0;
      break;
    default:
      assert(0);
    }
  }

  assert(intervals_size < intervals_capacity);

  lp_feasibility_set_t* result = lp_feasibility_set_new_from_intervals(intervals, intervals_size);

  // Construct the status
  if (all_s1) {
    *status = LP_FEASIBILITY_SET_INTERSECT_S1;
  } else if (all_s2) {
    *status = LP_FEASIBILITY_SET_INTERSECT_S2;
  } else if (result->size == 0) {
    *status = LP_FEASIBILITY_SET_EMPTY;
  } else {
    *status = LP_FEASIBILITY_SET_NEW;
  }

  for (i = intervals_size; i < intervals_capacity; ++ i) {
    lp_interval_destruct(intervals + i);
  }
  free(intervals);

  return result;
}

int lp_feasibility_set_empty(const lp_feasibility_set_t* set) {
  return set->size == 0;
}
