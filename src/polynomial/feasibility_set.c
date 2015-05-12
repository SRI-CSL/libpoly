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

#include "polynomial/polynomial.h"
#include "utils/debug_trace.h"

/**
 * Represents either an open interval or a point.
 */
typedef struct {

  enum {
    /** Interval is a point */
    FEASIBILITY_POINT,
    /** Interval is open (a, b) with a < b */
    FEASIBILITY_INTERVAL
  } type;

  /** The data itself */
  union {
    /** Value of the point */
    lp_value_t point_value;
    /** Value of the inteval */
    struct {
      lp_value_t lower_bound;
      lp_value_t upper_bound;
    } interval;
  };

} feasibility_interval_t;

/** Construct an interval (-inf, +inf) */
void feasibility_interval_construct(feasibility_interval_t* fi) {
  fi->type = FEASIBILITY_INTERVAL;
  lp_value_construct(&fi->interval.lower_bound, LP_VALUE_MINUS_INFINITY, 0);
  lp_value_construct(&fi->interval.upper_bound, LP_VALUE_PLUS_INFINITY, 0);
}

void feasibility_interval_destruct(feasibility_interval_t* fi) {
  switch (fi->type) {
  case FEASIBILITY_POINT:
    lp_value_destruct(&fi->point_value);
    break;
  case FEASIBILITY_INTERVAL:
    lp_value_destruct(&fi->interval.lower_bound);
    lp_value_destruct(&fi->interval.upper_bound);
  }
}

/** Returns true if the interval contains the value */
int feasibility_interval_contains(const feasibility_interval_t* interval, const lp_value_t* value) {
  switch (interval->type) {
  case FEASIBILITY_POINT:
    /** Just compare to the point */
    return (lp_value_cmp(&interval->point_value, value) == 0);
  case FEASIBILITY_INTERVAL: {
    /** Compare to lower and upper bounds */
    int cmp_lower = lp_value_cmp(&interval->interval.lower_bound, value);
    if (cmp_lower >= 0) {
      return 0;
    }
    int cmp_upper = lp_value_cmp(&interval->interval.upper_bound, value);
    if (cmp_upper <= 0) {
      return 0;
    }
    return 1;
  }
  }
  return 0;
}


int feasibility_interval_print(const feasibility_interval_t* fi, FILE* out) {
  int ret = 0;
  switch (fi->type) {
  case FEASIBILITY_POINT:
    ret += lp_value_print(&fi->point_value, out);
    break;
  case FEASIBILITY_INTERVAL:
    ret += fprintf(out, "(");
    ret += lp_value_print(&fi->interval.lower_bound, out);
    ret += fprintf(out, ", ");
    ret += lp_value_print(&fi->interval.upper_bound, out);
    ret += fprintf(out, ")");
  }
  return ret;
}

/**
 * Set of disjoint intervals representing an algebraic set, ordered from
 * left to right (-inf to +inf).
 */
struct lp_feasibility_set_struct {

  /** Number of intervals */
  size_t size;

  /** Capacity of the intervals table */
  size_t capacity;

  /** Vector feasibility intervals */
  feasibility_interval_t* intervals;

};

#define FEASIBILITY_SET_INITIAL_CAPACITY 10

static
void lp_feasibility_set_construct(lp_feasibility_set_t* s) {
  s->size = 0;
  s->capacity = FEASIBILITY_SET_INITIAL_CAPACITY;
  s->intervals = malloc(s->capacity * sizeof(feasibility_interval_t));
}

static
void lp_feasibility_set_destruct(lp_feasibility_set_t* s) {
  size_t i;
  for (i = 0; i < s->size; ++ i) {
    feasibility_interval_destruct(s->intervals + i);
  }
  free(s->intervals);
}

lp_feasibility_set_t* lp_feasibility_set_new() {
  lp_feasibility_set_t* result = malloc(sizeof(lp_feasibility_set_t));
  lp_feasibility_set_construct(result);
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
    ret += feasibility_interval_print(&set->intervals[i], out);
  }
  ret += fprintf(out, " }");
  return ret;
}

int lp_feasibility_set_contains(const lp_feasibility_set_t* set, const lp_value_t* value) {
  // TODO: binary search
  for (size_t i = 0; i < set->size; ++ i) {
    if (feasibility_interval_contains(&set->intervals[i], value)) {
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
