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

#include <variable_db.h>

#include "polynomial/output.h"
#include "utils/output.h"

#include <stdlib.h>
#include <string.h>

int monomial_print(const lp_polynomial_context_t* ctx, const lp_monomial_t* m, FILE* out) {
  int ret = 0;
  ret += integer_print(&m->a, out);
  ret += fprintf(out, " * ");
  size_t i = 0;
  for (i = 0; i < m->n; ++ i) {
    if (i) {
      ret += fprintf(out, "*");
    }
    ret += fprintf(out, "%s%s%zu", lp_variable_db_get_name(ctx->var_db, m->p[i].x), get_power_symbol(), m->p[i].d);
  }
  return ret;
}

int coefficient_print(const lp_polynomial_context_t* ctx, const coefficient_t* C, FILE* out) {
  int i, k = 0, ret = 0;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    ret += integer_print(&C->value.num, out);
    break;
  case COEFFICIENT_POLYNOMIAL: {
    // The polynomial
    const char* var_name = lp_variable_db_get_name(ctx->var_db, C->value.rec.x);
    for (i = SIZE(C) - 1; i >= 0; -- i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        switch (COEFF(C, i)->type) {
        case COEFFICIENT_POLYNOMIAL:

          if (k ++ ) {
            ret += fprintf(out, " + ");
          }
          ret += fprintf(out, "(");
          ret += coefficient_print(ctx, COEFF(C, i), out);
          ret += fprintf(out, ")");

          break;

        case COEFFICIENT_NUMERIC:

          if (integer_sgn(ctx->K, &COEFF(C, i)->value.num) > 0) {
            if (k ++) {
              ret += fprintf(out, " + ");
            }
            ret += integer_print(&COEFF(C, i)->value.num, out);
          } else {
            if (k ++) {
              ret += fprintf(out, " - ");
              lp_integer_t tmp;
              integer_construct_from_int(ctx->K, &tmp, 0);
              integer_neg(ctx->K, &tmp, &COEFF(C, i)->value.num);
              ret += integer_print(&tmp, out);
              integer_destruct(&tmp);
            } else {
              ret += integer_print(&COEFF(C, i)->value.num, out);
            }
          }

          break;
        }

        // Power
        if (i > 0) {
          if (i == 1) {
            ret += fprintf(out, "*%s", var_name);
          } else {
            ret += fprintf(out, "*%s%s%d", var_name, get_power_symbol(), i);
          }
        }
      }
    }
    break;
  }
  }
  return ret;
}

#if _XOPEN_SOURCE >= 700 || _POSIX_C_SOURCE >= 200809L
char* coefficient_to_string(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  coefficient_print(ctx, C, f);
  fclose(f);
  return str;
}
#endif
