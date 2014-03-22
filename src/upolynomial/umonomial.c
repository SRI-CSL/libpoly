/*
 * umonomial.c
 *
 *  Created on: Mar 22, 2014
 *      Author: dejan
 */

#include "upolynomial/umonomial.h"

#include <assert.h>

void umonomial_construct(int_ring K, umonomial_t* m, size_t degree, const integer_t* coefficient) {
  m->degree = degree;
  integer_ops.construct_copy(K, &m->coefficient, coefficient);
  assert(m->degree == 0 || integer_ops.sgn(Z, &m->coefficient));
}

void umonomial_construct_from_int(int_ring K, umonomial_t* m, size_t degree, long coefficient) {
  m->degree = degree;
  integer_ops.construct_from_int(K, &m->coefficient, coefficient);
  assert(m->degree == 0 || integer_ops.sgn(Z, &m->coefficient));
}

void umonomial_construct_copy(int_ring K, umonomial_t* m, const umonomial_t* from) {
  assert(from);
  umonomial_construct(K, m, from->degree, &from->coefficient);
  assert(m->degree == 0 || integer_ops.sgn(Z, &m->coefficient));
}

void umonomial_destruct(umonomial_t* m) {
  integer_ops.destruct(&m->coefficient);
}

int umonomial_print(const umonomial_t* m, FILE* out) {
  int len = 0;
  int sgn = integer_ops.sgn(Z, &m->coefficient);
  if (sgn < 0) len += fprintf(out, "(");
  len += integer_ops.print(&m->coefficient, out);
  if (m->degree) {
    if (m->degree == 1) {
      len += fprintf(out, "*x");
    } else {
      len += fprintf(out, "*x^%zu", m->degree);
    }
  }
  if (sgn < 0) len += fprintf(out, ")");
  return len;
}
