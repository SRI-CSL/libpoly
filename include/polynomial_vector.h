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

/** Allocate and construct a new vector */
lp_polynomial_vector_t* lp_polynomial_vector_new(const lp_polynomial_context_t* ctx);

/* Delete the vector */
void lp_polynomial_vector_delete(lp_polynomial_vector_t* v);

/** Add to back (makes a copy, should be in the context of the vector) */
void lp_polynomial_vector_push_back(lp_polynomial_vector_t* v, const lp_polynomial_t* p);

/** Reset the vector to 0 elements */
void lp_polynomial_vector_reset(lp_polynomial_vector_t* v);

/** Returns the size of the vector */
size_t lp_polynomial_vector_size(const lp_polynomial_vector_t* v);

/** Returns the polynomial at i (newly constructed each time) */
lp_polynomial_t* lp_polynomial_vector_at(const lp_polynomial_vector_t* v, size_t i);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
