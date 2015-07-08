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

#include <stdio.h>

/** Construct and attach an empty database */
lp_variable_db_t* lp_variable_db_new(void);

/** Attach the variable database */
void lp_variable_db_attach(lp_variable_db_t* var_db);

/** Detach the variable database */
void lp_variable_db_detach(lp_variable_db_t* var_db);

/** Add a new variable to the database and return its id */
lp_variable_t lp_variable_db_new_variable(lp_variable_db_t* var_db, const char* name);

/** Add a new variable to database with a given id */
void lp_variable_db_add_variable(lp_variable_db_t* var_db, lp_variable_t var, const char* name);

/** Print the variable database */
int lp_variable_db_print(const lp_variable_db_t* var_db, FILE* out);

/** Get the name of the variable */
const char* lp_variable_db_get_name(const lp_variable_db_t* var_db, lp_variable_t var);
