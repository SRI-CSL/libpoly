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

#include <upolynomial_factors.h>

#include "upolynomial/upolynomial.h"
#include "upolynomial/factors.h"
#include "upolynomial/output.h"

lp_upolynomial_factors_t* lp_upolynomial_factors_construct(void) {
  lp_upolynomial_factors_t* f = malloc(sizeof(lp_upolynomial_factors_t));
  integer_construct_from_int(lp_Z, &f->constant, 1);
  f->size = 0;
  f->capacity = 10;
  f->factors = calloc(f->capacity, sizeof(lp_upolynomial_t*));
  f->multiplicities = calloc(f->capacity, sizeof(size_t));
  return f;
}

void lp_upolynomial_factors_swap(lp_upolynomial_factors_t* f1, lp_upolynomial_factors_t* f2) {
  lp_upolynomial_factors_t tmp = *f1;
  *f1 = *f2;
  *f2 = tmp;
}

void lp_upolynomial_factors_clear(lp_upolynomial_factors_t* f) {
  size_t i;
  integer_assign_int(lp_Z, &f->constant, 1);
  for (i = 0; i < f->size; ++i) {
    if (f->factors[i]) {
      lp_upolynomial_delete(f->factors[i]);
    }
    f->factors[i] = 0;
  }
  f->size = 0;
}

void lp_upolynomial_factors_destruct(lp_upolynomial_factors_t* f, int destruct_factors) {
  if (destruct_factors) {
    lp_upolynomial_factors_clear(f);
  }
  integer_destruct(&f->constant);
  free(f->factors);
  free(f->multiplicities);
  free(f);
}

size_t lp_upolynomial_factors_size(const lp_upolynomial_factors_t* f) {
  return f->size;
}

lp_upolynomial_t* lp_upolynomial_factors_get_factor(lp_upolynomial_factors_t* f, size_t i, size_t* d) {
  *d = f->multiplicities[i];
  return f->factors[i];
}

const lp_integer_t* lp_upolynomial_factors_get_constant(const lp_upolynomial_factors_t* f) {
  return &f->constant;
}

void lp_upolynomial_factors_add(lp_upolynomial_factors_t* f, lp_upolynomial_t* p, size_t d) {
  // assert(upolynomial_degree(p) > 0); (we reuse this as general sets)

  if (f->size == f->capacity) {
    f->capacity *= 2;
    f->factors = realloc(f->factors, f->capacity*sizeof(lp_upolynomial_t*));
    f->multiplicities = realloc(f->multiplicities, f->capacity*sizeof(size_t));
  }
  f->factors[f->size] = p;
  f->multiplicities[f->size] = d;
  f->size ++;
}

const lp_int_ring_t* lp_upolynomial_factors_ring(const lp_upolynomial_factors_t* f) {
  if (f->size == 0) {
    return lp_Z;
  } else {
    return f->factors[0]->K;
  }
}

void lp_upolynomial_factors_set_ring(lp_upolynomial_factors_t* f, const lp_int_ring_t* K) {
  size_t i;
  for (i = 0; i < f->size; ++ i) {
    lp_upolynomial_set_ring(f->factors[i], K);
  }
}
