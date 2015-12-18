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

#include "polynomial_vector.h"
#include "polynomial.h"

#include <stdlib.h>

#define DEFAULT_SIZE 10

lp_polynomial_vector_t* lp_polynomial_vector_new(const lp_polynomial_context_t* ctx) {
  lp_polynomial_vector_t* result = (lp_polynomial_vector_t*) malloc(sizeof(lp_polynomial_vector_t));
  lp_polynomial_vector_construct(result, ctx);
  return result;
}

void lp_polynomial_vector_delete(lp_polynomial_vector_t* v) {
  lp_polynomial_vector_destruct(v);
  free(v);
}

void lp_polynomial_vector_construct(lp_polynomial_vector_t* v, const lp_polynomial_context_t* ctx) {
  v->ctx = ctx;
  v->capacity = DEFAULT_SIZE;
  v->size = 0;
  v->data = (coefficient_t*) malloc(DEFAULT_SIZE*sizeof(coefficient_t));
  lp_polynomial_context_attach((lp_polynomial_context_t*)ctx);
}

void lp_polynomial_vector_destruct(lp_polynomial_vector_t* v) {
  lp_polynomial_vector_reset(v);
  free(v->data);
  lp_polynomial_context_detach((lp_polynomial_context_t*)v->ctx);
}

static inline
void lp_polynomial_vector_check_size_for_add(lp_polynomial_vector_t* v) {
  if (v->size == v->capacity) {
    v->capacity ++;
    v->capacity += v->capacity >> 1;
    v->data = (coefficient_t*) realloc(v->data, v->capacity*sizeof(coefficient_t*));
  }
}

void lp_polynomial_vector_push_back(lp_polynomial_vector_t* v, const lp_polynomial_t* p) {
  lp_polynomial_vector_check_size_for_add(v);
  coefficient_construct_copy(v->ctx, v->data + v->size, &p->data);
  v->size ++;
}

void lp_polynomial_vector_push_back_coeff(lp_polynomial_vector_t* v, const coefficient_t* C) {
  lp_polynomial_vector_check_size_for_add(v);
  coefficient_construct_copy(v->ctx, v->data + v->size, C);
  v->size ++;
}

void lp_polynomial_vector_reset(lp_polynomial_vector_t* v) {
  size_t i = 0;
  for (i = 0; i < v->size; ++ i) {
    coefficient_destruct(v->data + i);
  }
  v->size = 0;
}

size_t lp_polynomial_vector_size(const lp_polynomial_vector_t* v) {
  return v->size;
}

lp_polynomial_t* lp_polynomial_vector_at(const lp_polynomial_vector_t* v, size_t i) {
  return lp_polynomial_new_from_coefficient(v->ctx, v->data + i);
}
