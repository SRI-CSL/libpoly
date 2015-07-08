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

#include <polynomial.h>

#include "polynomial/coefficient.h"

struct lp_polynomial_struct {
  /** The actual polynomial representation (so we can use it as a coefficient) */
  coefficient_t data;
  /** Hash */
  size_t hash;
  /** Is this an external polynomial (needs checks on function entry) */
  char external;
  /** Context of the polynomial */
  const lp_polynomial_context_t* ctx;
};

/** Construct from coefficient */
void lp_polynomial_construct_from_coefficient(lp_polynomial_t* A, const lp_polynomial_context_t* ctx, const coefficient_t* from);
