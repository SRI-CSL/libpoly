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

#include <interval.h>
#include <upolynomial.h>
#include <algebraic_number.h>

#include "upolynomial/upolynomial_dense.h"

#include "number/dyadic_rational.h"


/**
 * Compute the Sturm sequence of f. A Sturm sequence is a sequence
 *
 *  S[0] = f
 *  S[1] = f'
 *
 *  a*S[i-2] = Q*S[i-1] + b*S[i]
 *
 * with a*b < 0.
 */
void upolynomial_compute_sturm_sequence(const lp_upolynomial_t* f, upolynomial_dense_t* S, size_t* size);

/**
 * Count the number of real roots that the polynomial f has in the given open
 * interval. The polynomial f should be square-free.
 */
int upolynomial_roots_count_sturm(const lp_upolynomial_t* f, const lp_rational_interval_t* interval);

/**
 * Isolate the root intervals of the polynomial f and construct the resulting
 * numbers into the given array (should be at least deg(f) size). The size
 * will be updated to the number of roots.
 */
void upolynomial_roots_isolate_sturm(const lp_upolynomial_t* f, lp_algebraic_number_t* roots, size_t* roots_size);

/**
* Finds the roots for a polynomial over Zp. Uses brute force or rabin root
* finding, depending on p.
*/
void upolynomial_roots_find_Zp(const lp_upolynomial_t* f, lp_integer_t** roots, size_t* roots_size);
