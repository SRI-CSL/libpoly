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

#include <value.h>

void lp_value_approx(const lp_value_t* v, lp_rational_interval_t* approx);

/**
 * Case v1 and v2 to the same type.
 *
 * @param v1, v2 the values to case
 * @param v1_tmp the value to use for allocating a new value
 * @param v1_to_use, v2_to_use the cast values to use
 */
int lp_value_to_same_type(const lp_value_t* v1, const lp_value_t* v2,
    lp_value_t* v1_new, lp_value_t* v2_new,
    const lp_value_t** v1_to_use, const lp_value_t** v2_to_use);
