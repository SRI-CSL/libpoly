/*
 * monomial.c
 *
 *  Created on: Feb 10, 2014
 *      Author: dejan
 */

#include <monomial.h>
#include <variable_db.h>
#include "polynomial/polynomial.h"

#include "number/integer.h"

#include <stdlib.h>
#include <assert.h>
#include <stdlib.h>

void lp_monomial_construct(const lp_polynomial_context_t* ctx, lp_monomial_t* m) {
  integer_construct_from_int(ctx->K, &m->a, 0);
  m->n = 0;
  m->capacity = 0;
  m->p = 0;
}

void lp_monomial_clear(const lp_polynomial_context_t* ctx, lp_monomial_t* m) {
  integer_assign_int(ctx->K, &m->a, 0);
  m->n = 0;
}

void lp_monomial_construct_copy(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_monomial_t* from, int sort) {
  integer_construct_copy(ctx->K, &m->a, &from->a);
  m->n = from->n;
  m->capacity = from->n;
  m->p = malloc(m->n*sizeof(power_t));
  size_t i, j;
  // Copy
  for (i = 0; i < m->n; ++ i) {
    m->p[i] = from->p[i];
  }
  // Sort top first (not too many variables, do naive)
  if (sort) {
    for (i = 0; i < m->n; ++i) {
      for (j = i + 1; j < m->n; ++j) {
        if (lp_variable_order_cmp(ctx->var_order, m->p[i].x, m->p[j].x) < 0) {
          power_t tmp = m->p[i];
          m->p[i] = m->p[j];
          m->p[j] = tmp;
        }
      }
    }
  }
}

void lp_monomial_destruct(lp_monomial_t* m) {
  integer_destruct(&m->a);
  if (m->p) {
    free(m->p);
    m->p = 0;
    m->n = 0;
    m->capacity = 0;
  }
}

void lp_monomial_set_coefficient(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_integer_t* a) {
  lp_integer_assign(ctx->K, &m->a, a);
}

void lp_monomial_assign(const lp_polynomial_context_t* ctx, lp_monomial_t* m, const lp_monomial_t* from, int sort) {
  if (m != from) {
    lp_monomial_destruct(m);
    lp_monomial_construct_copy(ctx, m, from, sort);
  }
}

void lp_monomial_push(lp_monomial_t* m, lp_variable_t x, size_t d) {
  if (m->n == m->capacity) {
    m->capacity += 5;
    m->p = realloc(m->p, m->capacity*sizeof(power_t));
  }
  assert(m->n < m->capacity);
  m->p[m->n].x = x;
  m->p[m->n].d = d;
  ++ (m->n);
}

void lp_monomial_pop(lp_monomial_t* m) {
  assert(m->n > 0);
  -- m->n;
}

#define SWAP(m1, m2) { lp_monomial_t tmp = m1; m1 = m2; m2 = tmp; }
#define MIN(x, y) (x < y ? x : y)

void lp_monomial_gcd(const lp_polynomial_context_t* ctx, lp_monomial_t* gcd, const lp_monomial_t* m1, const lp_monomial_t* m2) {

  assert(ctx->K == lp_Z);

  lp_monomial_t result;
  lp_monomial_construct(ctx, &result);

  // GCD of the coefficients
  integer_gcd_Z(&result.a, &m1->a, &m2->a);

  // GCD of the power
  size_t m1_i = 0, m2_i = 0;
  while (m1_i < m1->n && m2_i < m2->n) {
    // Only keep powers that are equal
    // Variables in the monomial go top to bottom
    int var_cmp = lp_variable_order_cmp(ctx->var_order, m1->p[m1_i].x, m2->p[m2_i].x);
    if (var_cmp == 0) {
      lp_variable_t x = m1->p[m1_i].x;
      size_t d = MIN(m1->p[m1_i].d, m2->p[m2_i].d);
      lp_monomial_push(&result, x, d);
      m1_i ++;
      m2_i ++;
    } else if (var_cmp > 0) {
      m1_i ++;
    } else if (var_cmp < 0) {
      m2_i ++;
    }
  }

  SWAP(result, *gcd);
  lp_monomial_destruct(&result);
}

int lp_monomial_print(const lp_polynomial_context_t* ctx, const lp_monomial_t* m, FILE* out) {
  size_t i;
  int ret = 0;
  ret += lp_integer_print(&m->a, out);
  ret += fprintf(out, "*");
  for (i = 0; i < m->n; ++ i) {
    ret += fprintf(out, "%s^%zu", lp_variable_db_get_name(ctx->var_db, m->p[i].x), m->p[i].d);
  }
  return ret;
}
