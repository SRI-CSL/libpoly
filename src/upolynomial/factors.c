/*
 * factors.c
 *
 *  Created on: Mar 25, 2014
 *      Author: dejan
 */

#include "upolynomial/upolynomial.h"
#include "upolynomial/factors.h"

upolynomial_factors_t* upolynomial_factors_construct(void) {
  upolynomial_factors_t* f = malloc(sizeof(upolynomial_factors_t));
  integer_construct_from_int(Z, &f->constant, 1);
  f->size = 0;
  f->capacity = 10;
  f->factors = calloc(f->capacity, sizeof(upolynomial_t*));
  f->multiplicities = calloc(f->capacity, sizeof(size_t));
  return f;
}

void upolynomial_factors_swap(upolynomial_factors_t* f1, upolynomial_factors_t* f2) {
  upolynomial_factors_t tmp = *f1;
  *f1 = *f2;
  *f2 = tmp;
}

void upolynomial_factors_clear(upolynomial_factors_t* f) {
  size_t i;
  integer_assign_int(Z, &f->constant, 1);
  for (i = 0; i < f->size; ++i) {
    if (f->factors[i]) {
      upolynomial_destruct(f->factors[i]);
    }
  }
  f->size = 0;
}

void upolynomial_factors_destruct(upolynomial_factors_t* f, int destruct_factors) {
  if (destruct_factors) {
    upolynomial_factors_clear(f);
  }
  integer_destruct(&f->constant);
  free(f->factors);
  free(f->multiplicities);
  free(f);
}

size_t upolynomial_factors_size(const upolynomial_factors_t* f) {
  return f->size;
}

upolynomial_t* upolynomial_factors_get_factor(upolynomial_factors_t* f, size_t i, size_t* d) {
  *d = f->multiplicities[i];
  return f->factors[i];
}

const integer_t* upolynomial_factors_get_constant(const upolynomial_factors_t* f) {
  return &f->constant;
}

void upolynomial_factors_add(upolynomial_factors_t* f, upolynomial_t* p, size_t d) {
  // assert(upolynomial_degree(p) > 0); (we reuse this as general sets)

  if (f->size == f->capacity) {
    f->capacity *= 2;
    f->factors = realloc(f->factors, f->capacity*sizeof(upolynomial_t*));
    f->multiplicities = realloc(f->multiplicities, f->capacity*sizeof(size_t));
  }
  f->factors[f->size] = p;
  f->multiplicities[f->size] = d;
  f->size ++;
}

int upolynomial_factors_print(const upolynomial_factors_t* f, FILE* out) {
  int len = 0;
  len += integer_print(&f->constant, out);
  size_t i;
  for (i = 0; i < f->size; ++ i) {
    len += fprintf(out, " * ");
    len += fprintf(out, "[");
    len += upolynomial_print(f->factors[i], out);
    len += fprintf(out, "]^%zu", f->multiplicities[i]);
  }
  return len;
}

int_ring upolynomial_factors_ring(const upolynomial_factors_t* f) {
  if (f->size == 0) {
    return Z;
  } else {
    return f->factors[0]->K;
  }
}

void upolynomial_factors_set_ring(upolynomial_factors_t* f, int_ring K) {
  size_t i;
  for (i = 0; i < f->size; ++ i) {
    upolynomial_set_ring(f->factors[i], K);
  }
}

const upolynomial_factors_ops_t upolynomial_factors_ops = {
    upolynomial_factors_construct,
    upolynomial_factors_destruct,
    upolynomial_factors_clear,
    upolynomial_factors_swap,
    upolynomial_factors_size,
    upolynomial_factors_get_factor,
    upolynomial_factors_get_constant,
    upolynomial_factors_add,
    upolynomial_factors_print,
    upolynomial_factors_ring,
    upolynomial_factors_set_ring
};
