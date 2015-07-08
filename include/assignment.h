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
#include "value.h"

#include <stdio.h>

struct lp_assignment_struct {
  /** Size of the map */
  size_t size;
  /** The values */
  lp_value_t* values;
  /** The variable database */
  const lp_variable_db_t* var_db;
};

/** Construct an empty assignment */
void lp_assignment_construct(lp_assignment_t* m, const lp_variable_db_t* var_db);

/** Construct an empty assignment */
lp_assignment_t* lp_assignment_new(const lp_variable_db_t* var_db);

/** Destruct the assignment */
void lp_assignment_destruct(lp_assignment_t* m);

/** Destruct and free the assignment */
void lp_assignment_delete(lp_assignment_t* m);

/** Print the model */
int lp_assignment_print(const lp_assignment_t* m, FILE* out);

#if _XOPEN_SOURCE >= 700 || _POSIX_C_SOURCE >= 200809L
/** Get the string representation of the model */
char* lp_assignment_to_string(const lp_assignment_t* m);
#endif

/**
 * Set the value of a variable (value is copied over). If value is 0 (pointer)
 * the value is unset.
 */
void lp_assignment_set_value(lp_assignment_t* m, lp_variable_t x, const lp_value_t* value);

/** Get the value of a variable */
const lp_value_t* lp_assignment_get_value(const lp_assignment_t* m, lp_variable_t x);

/** Get an approximate value of the variable */
void lp_assignment_get_value_approx(const lp_assignment_t* m, lp_variable_t x, lp_rational_interval_t* approx);

/** Get the sign of the polynomial in the model */
int lp_assignment_sgn(const lp_assignment_t* m, const lp_polynomial_t* A);
