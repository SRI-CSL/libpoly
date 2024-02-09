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

#include <poly.h>

struct lp_upolynomial_vector_struct {
  /** Number of elements in the vector */
  size_t size;
  /** Size of the array */
  size_t capacity;
  /** The elements */
  lp_upolynomial_t**data;
};

typedef struct lp_upolynomial_vector_struct lp_upolynomial_vector_t;

/**
 * Construct a new upolynomial vector
 */
lp_upolynomial_vector_t* lp_upolynomial_vector_construct(void);

/**
 * Swaps two polynomial vectors
 */
void lp_upolynomial_vector_swap(lp_upolynomial_vector_t *v1, lp_upolynomial_vector_t *v2);

/**
 * Clears the upolynomial vector
 */
void lp_upolynomial_vector_clear(lp_upolynomial_vector_t *v);

/**
 * Deletes the vector.
 */
void lp_upolynomial_vector_delete(lp_upolynomial_vector_t *v);

/**
 * Returns the size of the vector.
 */
size_t lp_upolynomial_vector_size(const lp_upolynomial_vector_t *v);

/**
 * Returns a copy of the polynomial at position i
 */
lp_upolynomial_t* lp_upolynomial_vector_at(lp_upolynomial_vector_t *v, size_t i);

/**
 * Moves a polynomial to the vector. p must not be used afterwards.
 */
void lp_upolynomial_vector_move_back(lp_upolynomial_vector_t *v, lp_upolynomial_t* p);

/**
 * Makes a copy of p and adds it to the vector.
 */
void lp_upolynomial_vector_push_back(lp_upolynomial_vector_t *v, const lp_upolynomial_t* p);

/**
 * Removes and returns the last element of the vector.
 */
lp_upolynomial_t* lp_upolynomial_vector_pop(lp_upolynomial_vector_t *v);
