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

#include <stdlib.h>

#include "upolynomial_vector.h"

#include <upolynomial.h>
#include "upolynomial/upolynomial.h"

lp_upolynomial_vector_t *lp_upolynomial_vector_construct(void) {
  lp_upolynomial_vector_t*v = malloc(sizeof(lp_upolynomial_vector_t));
  v->size = 0;
  v->capacity = 10;
  v->data = malloc(v->capacity * sizeof(lp_upolynomial_t*));
  return v;
}

void lp_upolynomial_vector_swap(lp_upolynomial_vector_t *v1, lp_upolynomial_vector_t *v2) {
  lp_upolynomial_vector_t tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

void lp_upolynomial_vector_clear(lp_upolynomial_vector_t *v) {
  for (size_t i = 0; i < v->size; ++i) {
    lp_upolynomial_delete(v->data[i]);
  }
  v->size = 0;
}

void lp_upolynomial_vector_delete(lp_upolynomial_vector_t *v) {
  lp_upolynomial_vector_clear(v);
  free(v->data);
  free(v);
}

size_t lp_upolynomial_vector_size(const lp_upolynomial_vector_t *v) {
  return v->size;
}

lp_upolynomial_t *lp_upolynomial_vector_at(lp_upolynomial_vector_t *v, size_t i) {
  return lp_upolynomial_construct_copy(v->data[i]);
}

void lp_upolynomial_vector_move_back(lp_upolynomial_vector_t *v, lp_upolynomial_t *p) {
  if (v->size == v->capacity) {
    v->capacity *= 2;
    v->data = realloc(v->data, v->capacity * sizeof(lp_upolynomial_t*));
  }
  v->data[v->size] = p;
  v->size ++;
}

void lp_upolynomial_vector_push_back(lp_upolynomial_vector_t *v, const lp_upolynomial_t *p) {
  lp_upolynomial_vector_move_back(v, lp_upolynomial_construct_copy(p));
}

lp_upolynomial_t *lp_upolynomial_vector_pop(lp_upolynomial_vector_t *v) {
  if (v->size == 0) {
    return NULL;
  } else {
    v->size --;
    return v->data[v->size];
  }
}
