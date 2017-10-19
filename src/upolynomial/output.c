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

#include "upolynomial/output.h"
#include "upolynomial/factors.h"
#include "upolynomial/upolynomial.h"

#include "utils/output.h"

#include <stdio.h>

int umonomial_print(const ulp_monomial_t* m, FILE* out) {
  int len = 0;
  int sgn = integer_sgn(lp_Z, &m->coefficient);
  if (sgn < 0) len += fprintf(out, "(");
  len += integer_print(&m->coefficient, out);
  if (m->degree) {
    if (m->degree == 1) {
      len += fprintf(out, "*%s", get_upolynomial_var_symbol());
    } else {
      len += fprintf(out, "*x%s%zu", get_power_symbol(), m->degree);
    }
  }
  if (sgn < 0) len += fprintf(out, ")");
  return len;
}

int upolynomial_dense_print(const upolynomial_dense_t* p_d, FILE* file) {
  int len = 0;
  int k = p_d->size - 1;
  for (; k >= 0; --k) {
    int sgn = integer_sgn(lp_Z, p_d->coefficients + k);
    if (sgn) {
      if (sgn > 0) {
        fprintf(file, "+");
      }
      len += integer_print(p_d->coefficients + k, file);
      len += fprintf(file, "*x%s%d ", get_power_symbol(), k);
    }
  }
  return len;
}

int lp_upolynomial_print(const lp_upolynomial_t* p, FILE* out) {
  assert(p);

  int len = 0;
  unsigned i;
  for (i = 0; i < p->size; ++ i) {
    if (i) {
      len += fprintf(out, " + ");
    }
    len += umonomial_print(&p->monomials[p->size-i-1], out);
  }
  if (p->K) {
    len += fprintf(out, " [");
    len += lp_int_ring_print(p->K, out);
    len += fprintf(out, "]");
  }
  return len;
}

char* lp_upolynomial_to_string(const lp_upolynomial_t* p) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_upolynomial_print(p, f);
  fclose(f);
  return str;
}

int lp_upolynomial_factors_print(const lp_upolynomial_factors_t* f, FILE* out) {
  int len = 0;
  len += integer_print(&f->constant, out);
  size_t i;
  for (i = 0; i < f->size; ++ i) {
    len += fprintf(out, " * ");
    len += fprintf(out, "[");
    len += lp_upolynomial_print(f->factors[i], out);
    len += fprintf(out, "]^%zu", f->multiplicities[i]);
  }
  return len;
}
