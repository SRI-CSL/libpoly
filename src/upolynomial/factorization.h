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

#include <stdio.h>

#include <upolynomial.h>
#include <upolynomial_factors.h>

/**
 * Factors the given polynomial into a distinct degree factorization.
 * Polynomial f should in Z_p, square-free, and monic.
 *
 * The output of the function is a factorization where each factor is associated
 * with a distinct degree of each of its sub-factors.
 */
lp_upolynomial_factors_t* upolynomial_factor_distinct_degree(const lp_upolynomial_t* f);

/**
 * Factors the given polynomial using the algorithm of Berlekamp. The algorithm
 * assumes that p is in a ring Z_p for some prime p.
 */
lp_upolynomial_factors_t* upolynomial_factor_Zp(const lp_upolynomial_t* f);

/**
 * Factors the given polynomial using Hansel lifting.
 */
lp_upolynomial_factors_t* upolynomial_factor_Z(const lp_upolynomial_t* f);

/**
 * Factors a primitive polynomial into its square-free factors. If x is a factor,
 * it is factored into a separate factor.
 */
lp_upolynomial_factors_t* upolynomial_factor_square_free_primitive(const lp_upolynomial_t* f);
