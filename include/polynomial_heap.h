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

typedef int (*lp_polynomial_heap_compare_f)(const lp_polynomial_t* A, const lp_polynomial_t* B);

struct lp_polynomial_heap_struct {
    /** The data */
    lp_polynomial_t** data;
    /** Size of the data */
    size_t data_size;
    /** Number of set elements */
    size_t size;
    /** The compare function */
    lp_polynomial_heap_compare_f cmp;
};

/** Construct a new heap */
void lp_polynomial_heap_construct(lp_polynomial_heap_t* heap, lp_polynomial_heap_compare_f cmp);

/** Destruct the heap */
void lp_polynomial_heap_destruct(lp_polynomial_heap_t* heap);

/** Returns true if empty */
int lp_polynomial_heap_is_empty(lp_polynomial_heap_t* heap);

/** Add polynomial p to heap.  */
void lp_polynomial_heap_push(lp_polynomial_heap_t* heap, const lp_polynomial_t* p);

/** Removes and returns the top element of the heap, returns NULL if the heap is empty  */
lp_polynomial_t* lp_polynomial_heap_pop(lp_polynomial_heap_t* heap);

/** Clear the heap. */
void lp_polynomial_heap_clear(lp_polynomial_heap_t* heap);

/** Prints the heap */
void lp_polynomial_heap_print(lp_polynomial_heap_t* heap, FILE *);

#ifdef __cplusplus
} /* close extern "C" { */
#endif
