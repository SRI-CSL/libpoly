/**
 * Copyright 2024, SRI International.
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
#include "integer.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Set of disjoint intervals representing an algebraic set, ordered from
 * left to right (-inf to +inf).
 */
struct lp_feasibility_set_int_struct {
  /** The ring in which the feasibility set lives. */
  lp_int_ring_t *K;

  /** If true, set represents K \setminus elements */
  bool inverted;

  /** Number of elements */
  size_t size;

  /** Vector feasibility elements */
  lp_integer_t* elements;
};

/**
 * Create a new feasibility set K.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_new_full(lp_int_ring_t *K);

/**
 * Create a new feasibility set {}.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_new_empty(lp_int_ring_t *K);

/**
 * Construct a copy.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_new_copy(const lp_feasibility_set_int_t* set);

/**
 * Construct from integers.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_new_from_integer(lp_int_ring_t *K, const lp_integer_t* integers, size_t integers_size, bool inverted);

/**
 * Delete the given feasibility set.
 */
void lp_feasibility_set_int_delete(lp_feasibility_set_int_t* set);

/**
 * Assignment.
 */
void lp_feasibility_set_int_assign(lp_feasibility_set_int_t* set, const lp_feasibility_set_int_t* from);

/**
 * Swap.
 */
void lp_feasibility_set_int_swap(lp_feasibility_set_int_t* s1, lp_feasibility_set_int_t* s2);

/**
 * Check if the given set is empty.
 */
int lp_feasibility_set_int_is_empty(const lp_feasibility_set_int_t* set);

/**
 * Check if the given set is full.
 */
int lp_feasibility_set_int_is_full(const lp_feasibility_set_int_t* set);

/**
 * Check if the set is a point {a}.
 */
int lp_feasibility_set_int_is_point(const lp_feasibility_set_int_t* set);

/**
 * assigns the size of the set to out
 */
void lp_feasibility_set_int_size(const lp_feasibility_set_int_t *set, lp_integer_t *out);

/**
 * Check if the given value belongs to the set.
 */
int lp_feasibility_set_int_contains(const lp_feasibility_set_int_t* set, const lp_integer_t* value);

/**
 * Pick a value from the feasible set (must be non-empty).
 */
void lp_feasibility_set_int_pick_value(const lp_feasibility_set_int_t* set, lp_integer_t* value);

/**
 * Get intersection of the two sets.
 * s1 and s2 must be over the same ring K.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_intersect(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2);

/**
 * Get union of the two sets.
 * s1 and s2 must be over the same ring K.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_union(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2);

typedef enum {
  LP_FEASIBILITY_SET_INT_S1,
  LP_FEASIBILITY_SET_INT_S2,
  LP_FEASIBILITY_SET_INT_NEW,
  LP_FEASIBILITY_SET_INT_EMPTY
} lp_feasibility_set_int_status_t;

/**
 * Get intersection of the two sets, returns the status in the given variable.
 * The set s1 is given precedence so LP_FEASIBILITY_SET_S2 is the
 * status only if the intersect is not s1.
 * s1 and s2 must be over the same ring K.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_intersect_with_status(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2, lp_feasibility_set_int_status_t * status);

/**
 * Get union of the two sets, returns the status in the given variable.
 * The set s1 is given precedence so LP_FEASIBILITY_SET_S2 is the
 * status only if the union is not s1.
 * s1 and s2 must be over the same ring K.
 */
lp_feasibility_set_int_t* lp_feasibility_set_int_union_with_status(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2, lp_feasibility_set_int_status_t* status);

/**
 * Add one set to another, i.e. s = s \cup from.
 */
void lp_feasibility_set_int_add(lp_feasibility_set_int_t* s, const lp_feasibility_set_int_t* from);

/**
 * Print the set.
 */
int lp_feasibility_set_int_print(const lp_feasibility_set_int_t* set, FILE* out);

/**
 * Return the string representation of the set.
 */
char* lp_feasibility_set_int_to_string(const lp_feasibility_set_int_t* set);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
