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

/** Compute the greatest common divisor gcd(C1, C2). */
void coefficient_gcd(const lp_polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2);

/**
 * Compute the primitive part pp(A(divide with content).
 */
void coefficient_pp(const lp_polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C);

/**
 * Compute the content cont(A(gcd of coefficients). Content is computed so that
 * A/cont(A) has a positive leading coefficient.
 */
void coefficient_cont(const lp_polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C);

/** Comput pp and cont at the same time */
void coefficient_pp_cont(const lp_polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C);

/** Compute the least common multiple. */
void coefficient_lcm(const lp_polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2);

/** Compute the assumptions that ensure gcd is same degree, and put the assumptions into the assumptions vector. */
lp_polynomial_vector_t* coefficient_mgcd(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2, const lp_assignment_t* m);
