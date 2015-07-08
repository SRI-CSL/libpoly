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

#include "polynomial/coefficient.h"
#include "utils/output.h"

/** Print the monomial */
int monomial_print(const lp_polynomial_context_t* ctx, const lp_monomial_t* m, FILE* out);

/** Prints the coefficient to the given stream. */
int coefficient_print(const lp_polynomial_context_t* ctx, const coefficient_t* C, FILE* out);

/** Returns the string representation of the coefficient. */
char* coefficient_to_string(const lp_polynomial_context_t* ctx, const coefficient_t* C);
