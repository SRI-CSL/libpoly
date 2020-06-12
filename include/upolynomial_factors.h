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
#include "integer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Construct the factors with no factors, and constant 1. */
lp_upolynomial_factors_t* lp_upolynomial_factors_construct(void);

/**
 * Free the memory and the factor polynomials. If destruct_factors is true
 * then the individual factors are also destructed (should be, unless you
 * copied the factors somewhere else).
 */
void lp_upolynomial_factors_destruct(lp_upolynomial_factors_t* f, int destruct_factors);

/** Clear the factor polynomials */
void lp_upolynomial_factors_clear(lp_upolynomial_factors_t* f);

/** Swap the two factorizations */
void lp_upolynomial_factors_swap(lp_upolynomial_factors_t* f1, lp_upolynomial_factors_t* f2);

/** Get the number of factors */
size_t lp_upolynomial_factors_size(const lp_upolynomial_factors_t* f);

/** Get a factor with the given index i < size() */
lp_upolynomial_t* lp_upolynomial_factors_get_factor(lp_upolynomial_factors_t* f, size_t i, size_t* multiplicity);

/** Returns the constant of the factorization */
const lp_integer_t* lp_upolynomial_factors_get_constant(const lp_upolynomial_factors_t* f);

/** Add a factor with the given degree */
void lp_upolynomial_factors_add(lp_upolynomial_factors_t* f, lp_upolynomial_t* p, size_t d);

/** Print the factors */
int lp_upolynomial_factors_print(const lp_upolynomial_factors_t* f, FILE* out);

/** Get the ring */
const lp_int_ring_t* lp_upolynomial_factors_ring(const lp_upolynomial_factors_t* f);

/**
 * Set the ring of all polynomials to K. This is only possible if K is
 * "larger" than the existing ring.
 */
void lp_upolynomial_factors_set_ring(lp_upolynomial_factors_t* f, const lp_int_ring_t* K);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
