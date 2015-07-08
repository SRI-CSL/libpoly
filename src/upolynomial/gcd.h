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

#include <upolynomial.h>

/**
 * Compute the GCD using Euclid's algorithm. Degree of A should be >= than the
 * degree of B. This one only works with exact division, i.e. in Z_p rings.
 *
 * If U and V are not 0, then the extended gcd is computed, i.e. we solve
 *
 *   U*A + V*B = gcd(A, B) = D
 *
 * with deg(U) < deg(B) - deg(D) and deg(V) < deg(A) - deg(D).
 */
lp_upolynomial_t* upolynomial_gcd_euclid(const lp_upolynomial_t* A, const lp_upolynomial_t* B, lp_upolynomial_t** U, lp_upolynomial_t** V);

/**
 * Compute the GCD using the Subresultant algorithm. Degree of A should be >= than the
 * degree of B.
 */
lp_upolynomial_t* upolynomial_gcd_subresultant(const lp_upolynomial_t* A, const lp_upolynomial_t* B);


/**
 * Heuristic GCD computation. Pick a large value v, evaluate A(v), B(v), and
 * compute g = gcd(A(v), B(v)). Then reconstruct gcd from g. The method returns
 * 0 if the gcd compuation failed (heuristic!). Only on Z[x] polynomials.
 *
 * @param tries how many attempts
 */
lp_upolynomial_t* upolynomial_gcd_heuristic(const lp_upolynomial_t* A, const lp_upolynomial_t* B, int attempts);
