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

#include <stdio.h>

struct coefficient_factors_struct {
  size_t size;
  size_t capacity;
  coefficient_t* factors;
  size_t* multiplicities;
};

typedef struct coefficient_factors_struct coefficient_factors_t;

void coefficient_factors_construct(coefficient_factors_t* factors);
void coefficient_factors_destruct(coefficient_factors_t* factors);
void coefficient_factors_add(const lp_polynomial_context_t* ctx, coefficient_factors_t* factors, const coefficient_t* C, size_t multiplicity);
int coefficient_factors_print(const lp_polynomial_context_t* ctx, const coefficient_factors_t* factors, FILE* out);

/** Factors the given coefficient into square-free factors. */
void coefficient_factor_square_free(const lp_polynomial_context_t* ctx, const coefficient_t* C, coefficient_factors_t* factors);

/** Factor the given coefficient into content-free factors */
void coefficient_factor_content_free(const lp_polynomial_context_t* ctx, const coefficient_t* C, coefficient_factors_t* factors);
