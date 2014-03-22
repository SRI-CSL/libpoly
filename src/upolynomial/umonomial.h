#pragma once

#include "number/integer.h"

/**
 * A monomial is a coefficient and the degree.
 */
typedef struct umonomial_struct {
  /** Degree of the monomial */
  size_t degree;
  /** Coefficient with the monomial */
  integer_t coefficient;
} umonomial_t;

void umonomial_construct(int_ring K, umonomial_t* m, size_t degree, const integer_t* coefficient);

void umonomial_construct_from_int(int_ring K, umonomial_t* m, size_t degree, long coefficient);

void umonomial_construct_copy(int_ring K, umonomial_t* m, const umonomial_t* from);

void umonomial_destruct(umonomial_t* m);

int umonomial_print(const umonomial_t* m, FILE* out);

