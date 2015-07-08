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
#include "upolynomial/upolynomial.h"

/**
 * Upper bound on the modulus of the roots of f.
 *
 * Given
 *
 *  f = a_n x^n + ... + a_0, with a_n != 0
 *
 * the bound is
 *
 *  B = 1 + max(a_0, ..., a_{n-1})/|a_n|
 */
void upolynomial_root_bound_cauchy(const lp_upolynomial_t* f, lp_integer_t* B);

/**
 * Computes the bound on the size of coefficients of any polynomial g with
 * deg(g) <= n that divides f.
 *
 * Given
 *
 *  f = a_n x^n + ... + a_0, with a_n != 0
 *
 * the bound is
 *
 *  B = 2^d*norm(f)
 */
void upolynomial_factor_bound_landau_mignotte(const lp_upolynomial_t* f, size_t n, lp_integer_t* B);
