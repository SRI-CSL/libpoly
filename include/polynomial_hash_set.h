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

#ifdef __cplusplus
extern "C" {
#endif

struct lp_polynomial_hash_set_struct {
  /** The data */
  lp_polynomial_t** data;
  /** Size of the data */
  size_t data_size;
  /** Number of set elements */
  size_t size;
  /** Treshold for resize */
  size_t resize_threshold;
  /** Has the set been closed */
  int closed;
};

/** Construct a new set */
void lp_polynomial_hash_set_construct(lp_polynomial_hash_set_t* set);

/** Destruct the set */
void lp_polynomial_hash_set_destruct(lp_polynomial_hash_set_t* set);

/** Returns true if empty */
int lp_polynomial_hash_set_is_empty(lp_polynomial_hash_set_t* set);

/** Check whether p is in set. The set must not be closed). */
int lp_polynomial_hash_set_contains(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p);

/** Add polynomial p to set. Returns true if p was added (not already in the set).  */
int lp_polynomial_hash_set_insert(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p);

/** Close the set: compact the data so that all elements get stored in data[0..size]. No addition after close! */
void lp_polynomial_hash_set_close(lp_polynomial_hash_set_t* set);

/** Clear the set. */
void lp_polynomial_hash_set_clear(lp_polynomial_hash_set_t* set);

/** Print the set. */
int lp_polynomial_hash_set_print(const lp_polynomial_hash_set_t* set, FILE* out);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
