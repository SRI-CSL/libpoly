#pragma once

#include "number/integer.h"

/**
 * A monomial is a coefficient and the degree.
 */
typedef struct umonomial_struct {
  /** Degree of the monomial */
  size_t degree;
  /** Coefficient with the monomial */
  lp_integer_t coefficient;
} ulp_monomial_t;

void umonomial_construct(lp_int_ring_t* K, ulp_monomial_t* m, size_t degree, const lp_integer_t* coefficient);

void umonomial_construct_from_int(lp_int_ring_t* K, ulp_monomial_t* m, size_t degree, long coefficient);

void umonomial_construct_copy(lp_int_ring_t* K, ulp_monomial_t* m, const ulp_monomial_t* from);

void umonomial_destruct(ulp_monomial_t* m);

