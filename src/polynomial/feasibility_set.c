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

#include <feasibility_set.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "poly.h"
#include "value.h"
#include "interval.h"

#include "polynomial/polynomial.h"
#include "utils/debug_trace.h"
#include "polynomial/feasibility_set.h"

static
void lp_feasibility_set_ensure_capacity(lp_feasibility_set_t* s, size_t capacity) {
  if (capacity && capacity > s->capacity) {
    s->capacity = capacity;
    s->intervals = realloc(s->intervals, s->capacity * sizeof(lp_interval_t));
  }
}

static
void lp_feasibility_set_construct(lp_feasibility_set_t* s, size_t size) {
  s->size = 0;
  s->capacity = 0;
  s->intervals = NULL;
  lp_feasibility_set_ensure_capacity(s, size);
}

lp_feasibility_set_t* lp_feasibility_set_new_empty(void) {
  return lp_feasibility_set_new_internal(0);
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

lp_feasibility_set_t* lp_feasibility_set_new_full(void) {
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

static
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

void lp_feasibility_set_construct_from_interval(lp_feasibility_set_t* set, const lp_interval_t* from) {
  lp_feasibility_set_construct(set, 1);
  lp_interval_construct_copy(set->intervals, from);
  set->size = 1;
}

lp_feasibility_set_t* lp_feasibility_set_new_from_interval(const lp_interval_t* I) {
  lp_feasibility_set_t* result = malloc(sizeof(lp_feasibility_set_t));
  lp_feasibility_set_construct_from_interval(result, I);
  return result;
}

void lp_feasibiliy_set_assign(lp_feasibility_set_t* set, const lp_feasibility_set_t* from) {
  if (set != from) {
    lp_feasibility_set_destruct(set);
    lp_feasibility_set_construct_copy(set, from);
  }
}

void lp_feasibility_set_swap(lp_feasibility_set_t* s1, lp_feasibility_set_t* s2) {
  lp_feasibility_set_t tmp = *s1;
  *s1 = *s2;
  *s2 = tmp;
}

int lp_feasibility_set_is_empty(const lp_feasibility_set_t* set) {
  return set->size == 0;
}

int lp_feasibility_set_is_full(const lp_feasibility_set_t* set) {
  if (set->size != 1) {
    return 0;
  }
  return lp_interval_get_lower_bound(set->intervals)->type == LP_VALUE_MINUS_INFINITY &&
      lp_interval_get_upper_bound(set->intervals)->type == LP_VALUE_PLUS_INFINITY;
}

int lp_feasibility_set_is_point(const lp_feasibility_set_t* set) {
  return set->size == 1 && set->intervals->is_point;
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
  struct u_memstream mem;
  char* str = 0;
  size_t size = 0;
  u_memstream_open(&mem, &str, &size);
  FILE* f = u_memstream_get(&mem);
  lp_feasibility_set_print(set, f);
  u_memstream_close(&mem);
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

int lp_feasibility_set_contains_int(const lp_feasibility_set_t* set) {
  for (size_t i = 0; i < set->size; ++ i) {
    if (lp_interval_contains_int(set->intervals + i)) {
      return 1;
    }
  }
  return 0;
}

// We get smallest integer < rational < algebraic
// If same we get one with largest interval size
static inline
int same_or_better_complexity(const lp_value_t* v1, int v1_interval_size, const lp_value_t* v2, int v2_interval_size) {

  // See if any is integer
  int v1_is_int = lp_value_is_integer(v1);
  int v2_is_int = lp_value_is_integer(v2);
  if (v1_is_int && !v2_is_int) {
    return 1;
  }
  if (v2_is_int && !v1_is_int) {
    return 0;
  }

  // See if any is rational
  int v1_is_rational = lp_value_is_rational(v1);
  int v2_is_rational = lp_value_is_rational(v2);
  if (v1_is_rational && !v2_is_rational) {
    return 1;
  }
  if (v2_is_rational && !v1_is_rational) {
    return 0;
  }

  // Same type, compare the intervals
  return v1_interval_size >= v2_interval_size;
}

void lp_feasibility_set_pick_value(const lp_feasibility_set_t* set, lp_value_t* value) {
  size_t i;

  assert(!lp_feasibility_set_is_empty(set));

  lp_interval_pick_value(set->intervals, value);
  int value_interval_size = lp_interval_size_approx(set->intervals);

  lp_value_t current;
  lp_value_construct_none(&current);
  for (i = 1; i < set->size; ++ i) {
    int current_interval_size = lp_interval_size_approx(set->intervals + i);
    lp_interval_pick_value(set->intervals + i, &current);
    if (!same_or_better_complexity(value, value_interval_size, &current, current_interval_size)) {
      lp_value_swap(value, &current);
      value_interval_size = current_interval_size;
    }
  }
  lp_value_destruct(&current);
}

void lp_feasibility_set_pick_first_value(const lp_feasibility_set_t* set, lp_value_t* value) {
  assert(!lp_feasibility_set_is_empty(set));
  lp_interval_pick_value(set->intervals, value);
}

lp_feasibility_set_t* lp_feasibility_set_intersect(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2) {
  lp_feasibility_set_intersect_status_t status;
  return lp_feasibility_set_intersect_with_status(s1, s2, &status);
}

lp_feasibility_set_t* lp_feasibility_set_intersect_with_status(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2, lp_feasibility_set_intersect_status_t* status) {

  // Corner cases
  if (s1->size == 0 || s2->size == 0) {
    *status = LP_FEASIBILITY_SET_EMPTY;
    return lp_feasibility_set_new_internal(0);
  }

  // Size of the result is at most max of the sizes
  size_t intervals_capacity = s1->size + s2->size;
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

  // Whoever didn't run out is not all in
  if (s1_i < s1->size) {
    all_s1 = 0;
  }
  if (s2_i < s2->size) {
    all_s2 = 0;
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

/**
 * Sort the intervals so that we can scan from left to right and keep adding
 * one of the upper bounds to the union. Basically we're sorting by lower bound.
 */
static
int interval_sort_for_union(const void *I1_void, const void* I2_void) {
  const lp_interval_t* I1 = (const lp_interval_t*) I1_void;
  const lp_interval_t* I2 = (const lp_interval_t*) I2_void;
  lp_interval_cmp_t cmp = lp_interval_cmp(I1, I2);

  switch (cmp) {
  /* I1: (  )
   * I2:      (   ) */
  case LP_INTERVAL_CMP_LT_NO_INTERSECT:
    // Scan I1 first, then I2
    return -1;
  /* I1: (   )
   * I2:   (   )    */
  case LP_INTERVAL_CMP_LT_WITH_INTERSECT:
    // Scan I1 first, then I2
    return -1;
  /* I1: (   )
   * I2: (     )    */
  case LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1:
    // I2 might have smaller lower bound, so we scan it first
    return 1;
  /* I1: (     ]
   * I2:   (   ]    */
  case LP_INTERVAL_CMP_LEQ_WITH_INTERSECT_I2:
    // Scan I1 fits, then I2
    return -1;
  /* I1: (   ]
   * I2: (   ]      */
  case LP_INTERVAL_CMP_EQ:
    // Doesn't matter
    return 0;
  /* I1:   (   ]
   * I2: (     ]    */
  case LP_INTERVAL_CMP_GEQ_WITH_INTERSECT_I1:
    // Scan I2 first
    return 1;
  /* I1: (       )
   * I2: (    )     */
  case LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2:
    // Scan I1 first, it might have smaller lower bound
    return -1;
  /* I1:   (    )
   * I2: (    )     */
  case LP_INTERVAL_CMP_GT_WITH_INTERSECT:
    // Scan I2 first
    return 1;
  /* I1:      (   )
   * I2: (  )       */
  case LP_INTERVAL_CMP_GT_NO_INTERSECT:
    // Scan I2 first
    return 1;
  }

  return 1;
}

void lp_feasibility_set_add(lp_feasibility_set_t* s, const lp_feasibility_set_t* from) {

  if (lp_feasibility_set_is_empty(from)) {
    return;
  }

  if (lp_feasibility_set_is_full(s)) {
    return;
  }

  // Make sure there is space
  lp_feasibility_set_ensure_capacity(s, s->size + from->size);
  // Copy over the intervals
  size_t i;
  lp_interval_t* copy_to = s->intervals + s->size;
  for (i = 0; i < from->size; ++ i, copy_to ++) {
    lp_interval_construct_copy(copy_to, from->intervals + i);
  }
  s->size += from->size;
  // Add the new intervals
  qsort(s->intervals, s->size, sizeof(lp_interval_t), interval_sort_for_union);

  if (trace_is_enabled("feasibility_set")) {
    for (i = 0; i < s->size; ++ i) {
      lp_interval_print(s->intervals + i, trace_out); tracef("\n");
    }
  }

  // Now, normalize the intervals
  size_t keep;
  for (i = 1, keep = 1; i < s->size; ++ i) {
    // If new one has intersection with previous one, merge
    // (   )
    //    (  )
    // (     )
    // otherwise keep the interval

    // Compare and decide whether to merge
    const lp_interval_t* I1 = s->intervals + keep - 1;
    const lp_interval_t* I2 = s->intervals + i;
    lp_interval_cmp_t cmp = lp_interval_cmp(I1, I2);

    if (trace_is_enabled("feasibility_set")) {
      tracef("I1 = "); lp_interval_print(I1, trace_out); tracef("\n");
      tracef("I2 = "); lp_interval_print(I2, trace_out); tracef("\n");
    }

    int merge = 0;
    int ignore = 0;

    switch (cmp) {
    /* I1: (  )
     * I2:      (   ) */
    case LP_INTERVAL_CMP_LT_NO_INTERSECT:
      // Check if the edges are the same
      if (lp_value_cmp(lp_interval_get_upper_bound(I1), lp_interval_get_lower_bound(I2)) == 0 && (!I1->b_open || !I2->a_open)) {
        merge = 1;
      }
      break;
    /* I1: (   )
     * I2:   (   )    */
    case LP_INTERVAL_CMP_LT_WITH_INTERSECT:
      merge = 1;
      break;
    /* I1: (   )
     * I2: (     )    */
    case LP_INTERVAL_CMP_LT_WITH_INTERSECT_I1:
      // I2 should have been first
      assert(0);
      break;
    /* I1: (     ]
     * I2:   (   ]    */
    case LP_INTERVAL_CMP_LEQ_WITH_INTERSECT_I2:
      merge = 1;
      break;
    /* I1: (   ]
     * I2: (   ]      */
    case LP_INTERVAL_CMP_EQ:
      merge = 1;
      break;
    /* I1:   (   ]
     * I2: (     ]    */
    case LP_INTERVAL_CMP_GEQ_WITH_INTERSECT_I1:
      merge = 1;
      break;
    /* I1: (       )
     * I2: (    )     */
    case LP_INTERVAL_CMP_GT_WITH_INTERSECT_I2:
      // Just keep I1
      ignore = 1;
      break;
    /* I1:   (    )
     * I2: (    )     */
    case LP_INTERVAL_CMP_GT_WITH_INTERSECT:
      // I2 should have been first
      assert(0);
      break;
    /* I1:      (   )
     * I2: (  )       */
    case LP_INTERVAL_CMP_GT_NO_INTERSECT:
      // I2 should have been first
      assert(0);
      break;
    }

    // Merge if asked
    if (merge) {
      // Just use the endpoint
      const lp_value_t* new_b = lp_interval_get_upper_bound(s->intervals + i);
      lp_interval_set_b(s->intervals + keep - 1, new_b, s->intervals[i].b_open);
    } else if (!ignore) {
      // We just keep the new one
      if (i != keep) {
        lp_interval_swap(s->intervals + i, s->intervals + keep);
      }
      keep ++;
    }
  }

  // Destroy the leftover ones
  for (i = keep; i < s->size; ++ i) {
    lp_interval_destruct(s->intervals + i);
  }
  s->size = keep;

  if (trace_is_enabled("feasibility_set")) {
    for (i = 0; i < s->size; ++ i) {
      lp_interval_print(s->intervals + i, trace_out); tracef("\n");
    }
  }
}

void lp_feasibility_set_to_interval(const lp_feasibility_set_t* set, lp_interval_t* result) {
  assert(set);
  assert(set->size > 0);
  const lp_interval_t* first = set->intervals;
  const lp_value_t* a = &first->a;
  const lp_interval_t* last = set->intervals + set->size - 1;
  const lp_value_t* b = last->is_point ? &last->a : &last->b;
  lp_interval_t tmp;
  lp_interval_construct(&tmp, a, first->a_open, b, last->b_open);
  lp_interval_swap(result, &tmp);
  lp_interval_destruct(&tmp);
}
