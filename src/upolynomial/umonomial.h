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

#include "number/integer.h"

/**
 * A monomial is a coefficient and the degree.
 */
typedef struct umonomial_struct {
  /** Degree of the monomial */
  size_t degree;
  /** Coefficient with the monomial */
  lp_integer_t coefficient;
} ulp_monomial_t;

void umonomial_construct(const lp_int_ring_t* K, ulp_monomial_t* m, size_t degree, const lp_integer_t* coefficient);

void umonomial_construct_from_int(const lp_int_ring_t* K, ulp_monomial_t* m, size_t degree, long coefficient);

void umonomial_construct_copy(const lp_int_ring_t* K, ulp_monomial_t* m, const ulp_monomial_t* from);

void umonomial_destruct(ulp_monomial_t* m);
