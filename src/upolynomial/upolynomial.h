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

#include <integer.h>
#include "upolynomial/umonomial.h"

/**
 * A polynomial is the ring, number of monomials and the monomials.
 */
struct lp_upolynomial_struct {
  /** The ring of coefficients */
  lp_int_ring_t* K;
  /** The number of monomials */
  size_t size;
  /** The monomials */
  ulp_monomial_t monomials[];
};
