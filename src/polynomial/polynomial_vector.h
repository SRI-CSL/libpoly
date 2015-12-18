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

#include <polynomial_vector.h>

#include "coefficient.h"

struct lp_polynomial_vector_struct {
  /** The context */
  const lp_polynomial_context_t* ctx;
  /** Capacity of the vector */
  size_t capacity;
  /** Size of the vector */
  size_t size;
  /** The polynomials */
  coefficient_t* data;
};

/** Construct the vector */
void lp_polynomial_vector_construct(lp_polynomial_vector_t* v, const lp_polynomial_context_t* ctx);

/* Destruct the vector */
void lp_polynomial_vector_destruct(lp_polynomial_vector_t* v);

/** Add to back (makes a copy) */
void lp_polynomial_vector_push_back_coeff(lp_polynomial_vector_t* v, const coefficient_t* C);
