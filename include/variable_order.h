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

#ifdef __cplusplus
extern "C" {
#endif

/** Construct a new variable order. Attaches one. */
lp_variable_order_t* lp_variable_order_new(void);

/** Attach an object to this order (constructor should attach) */
void lp_variable_order_attach(lp_variable_order_t* var_order);

/** Detach an object from this order (frees if refcount = 0) */
void lp_variable_order_detach(lp_variable_order_t* var_order);

/**
 * Compare two variables:
 *  * by index in the push order
 *  * if the variable hasn't been pushed it's index is +inf
 *  * if indices are the same, compare by variable id
 *
 * In polynomial operations, the variables are order so that the highest
 * variable is the main variable.
  */
int lp_variable_order_cmp(const lp_variable_order_t* var_order, lp_variable_t x, lp_variable_t y);

/** returns the bigger of the two variables wrt. the order */
static inline
lp_variable_t lp_variable_order_max(const lp_variable_order_t* var_order, lp_variable_t x, lp_variable_t y) {
  if (x == lp_variable_null) return y;
  if (y == lp_variable_null) return x;
  return lp_variable_order_cmp(var_order, x, y) < 0 ? x : y;
}

/** returns the smaller of the two variables wrt. the order */
static inline
lp_variable_t lp_variable_order_min(const lp_variable_order_t* var_order, lp_variable_t x, lp_variable_t y) {
  if (x == lp_variable_null) return y;
  if (y == lp_variable_null) return x;
  return lp_variable_order_cmp(var_order, x, y) > 0 ? x : y;
}

/** Get the size of the order */
size_t lp_variable_order_size(const lp_variable_order_t* var_order);

/** Clear the order */
void lp_variable_order_clear(lp_variable_order_t* var_order);

/** Does the order have an opinion on x */
int lp_variable_order_contains(lp_variable_order_t* var_order, lp_variable_t x);

/** Push a variable to the list */
void lp_variable_order_push(lp_variable_order_t* var_order, lp_variable_t var);

/** Pop the last variable from the list */
void lp_variable_order_pop(lp_variable_order_t* var_order);

/** Reverse the order. This only affects the compare function, not push and pop */
void lp_variable_order_reverse(lp_variable_order_t* var_order);

/** Get the last variable from the list */
lp_variable_t lp_variable_order_top(const lp_variable_order_t* var_order);

/** Print the list of variables */
int lp_variable_order_print(const lp_variable_order_t* var_order, const lp_variable_db_t* var_db, FILE* out);

/** Return a string representation of the order */
char* lp_variable_order_to_string(const lp_variable_order_t* var_order, const lp_variable_db_t* var_db);

/** Return the variable list */
const lp_variable_list_t* lp_variable_order_get_list(const lp_variable_order_t* var_order);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
