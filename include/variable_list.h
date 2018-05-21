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

/**
 * A list of variable that keeps an index for each variable. We don't expect
 * too many variables in the database, so we keep a map as an array.
 */
struct lp_variable_list_struct {
  /** List of variables in the order */
  lp_variable_t *list;
  /** Size of the list */
  size_t list_size;
  /** Capacity of the list */
  size_t list_capacity;
  /** Map from variables to the index in the list (-1 if not in the order) */
  int* var_to_index_map;
  /** Size of the variable map */
  size_t var_to_index_map_capacity;
};

/** Construct a new variable order */
void lp_variable_list_construct(lp_variable_list_t* list);

/** Destruct the variable order */
void lp_variable_list_destruct(lp_variable_list_t* list);

/** Get the size of the list */
size_t lp_variable_list_size(const lp_variable_list_t* list);

/** Get the index of the variable in the list, or -1 if not there */
int lp_variable_list_index(const lp_variable_list_t* list, lp_variable_t x);

/** Copy the variables into the given vector */
void lp_variable_list_copy_into(const lp_variable_list_t* list, lp_variable_t* vars);

/** Push a variable to the list */
void lp_variable_list_push(lp_variable_list_t* list, lp_variable_t var);

/** Pop the last variable from the list */
void lp_variable_list_pop(lp_variable_list_t* list);

/** Get the last variable from the list */
lp_variable_t lp_variable_list_top(const lp_variable_list_t* list);

/** Returns 1 if the list contains the given variable */
int lp_variable_list_contains(const lp_variable_list_t* list, lp_variable_t x);

/** Removes a variable from the list, i.e. replace with lp_variable_null */
void lp_variable_list_remove(lp_variable_list_t* list, lp_variable_t x);

/** Order the list based on the given order */
void lp_variable_list_order(lp_variable_list_t* list, const lp_variable_order_t* order);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
