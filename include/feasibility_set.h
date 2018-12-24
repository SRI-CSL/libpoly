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
#include "sign_condition.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  lp_interval_t* intervals;

};

/**
 * Create a new feasibility set [-inf, inf].
 */
lp_feasibility_set_t* lp_feasibility_set_new_full(void);

/**
 * Create a new feasibility set {}.
 */
lp_feasibility_set_t* lp_feasibility_set_new_empty(void);

/**
 * Construct a copy.
 */
lp_feasibility_set_t* lp_feasibility_set_new_copy(const lp_feasibility_set_t* set);

/**
 * Delete the given feasibility set.
 */
void lp_feasibility_set_delete(lp_feasibility_set_t* set);

/**
 * Assignment.
 */
void lp_feasibiliy_set_assign(lp_feasibility_set_t* set, const lp_feasibility_set_t* from);

/**
 * Swap.
 */
void lp_feasibility_set_swap(lp_feasibility_set_t* s1, lp_feasibility_set_t* s2);

/**
 * Check if the given set is empty.
 */
int lp_feasibility_set_is_empty(const lp_feasibility_set_t* set);

/**
 * Check if the given set is full, i.e. (-inf, +inf).
 */
int lp_feasibility_set_is_full(const lp_feasibility_set_t* set);

/**
 * Check if the set is a point {a}.
 */
int lp_feasibility_set_is_point(const lp_feasibility_set_t* set);

/**
 * Check if the given value belongs to the set.
 */
int lp_feasibility_set_contains(const lp_feasibility_set_t* set, const lp_value_t* value);

/**
 * Pick a value from the feasible set (must be non-empty). If an integer value
 * is available it will be picked.
 */
void lp_feasibility_set_pick_value(const lp_feasibility_set_t* set, lp_value_t* v);

/**
 * Pick a value from the first interval towards -inf.
 */
void lp_feasibility_set_pick_first_value(const lp_feasibility_set_t* set, lp_value_t* v);

/**
 * Get intersection of the two sets, returns the status in the given variable.
 */
lp_feasibility_set_t* lp_feasibility_set_intersect(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2);

typedef enum {
  LP_FEASIBILITY_SET_INTERSECT_S1,
  LP_FEASIBILITY_SET_INTERSECT_S2,
  LP_FEASIBILITY_SET_NEW,
  LP_FEASIBILITY_SET_EMPTY
} lp_feasibility_set_intersect_status_t;

/**
 * Get intersection of the two sets, returns the status in the given variable.
 * The set s1 is given precedence so LP_FEASIBILITY_SET_INTERSECT_S2 is the
 * status only if the intersect is not s1.
 */
lp_feasibility_set_t* lp_feasibility_set_intersect_with_status(const lp_feasibility_set_t* s1, const lp_feasibility_set_t* s2, lp_feasibility_set_intersect_status_t* status);

/**
 * Add one set to another, i.e. s = s \cup from.
 */
void lp_feasibility_set_add(lp_feasibility_set_t* s, const lp_feasibility_set_t* from);


/**
 * Print the set.
 */
int lp_feasibility_set_print(const lp_feasibility_set_t* set, FILE* out);

/**
 * Return the string representation of the set.
 */
char* lp_feasibility_set_to_string(const lp_feasibility_set_t* set);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
