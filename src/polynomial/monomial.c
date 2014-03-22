/*
 * monomial.c
 *
 *  Created on: Feb 10, 2014
 *      Author: dejan
 */

#include "polynomial/monomial.h"
#include "polynomial/polynomial_internal.h"

#include "number/integer.h"

#include <malloc.h>
#include <assert.h>
#include <stdlib.h>

static
void monomial_construct(const polynomial_context_t* ctx, monomial_t* m) {
  integer_construct_from_int(ctx->K, &m->a, 0);
  m->n = 0;
  m->capacity = 0;
  m->p = 0;
}

static
void monomial_clear(const polynomial_context_t* ctx, monomial_t* m) {
  integer_assign_int(ctx->K, &m->a, 0);
  m->n = 0;
}

static
void monomial_construct_copy(const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort) {
  integer_construct_copy(ctx->K, &m->a, &from->a);
  m->n = from->n;
  m->capacity = from->n;
  m->p = malloc(m->n*sizeof(power_t));
  int i,j;
  // Copy
  for (i = 0; i < m->n; ++ i) {
    m->p[i] = from->p[i];
  }
  // Sort top first (not too many variables, do naive)
  if (sort) {
    for (i = 0; i < m->n; ++i) {
      for (j = i + 1; j < m->n; ++j) {
        if (ctx->var_order->ops->cmp(ctx->var_order, m->p[i].x, m->p[j].x) < 0) {
          power_t tmp = m->p[i];
          m->p[i] = m->p[j];
          m->p[j] = tmp;
        }
      }
    }
  }
}

static
void monomial_destruct(monomial_t* m) {
  integer_destruct(&m->a);
  if (m->p) {
    free(m->p);
    m->p = 0;
    m->n = 0;
    m->capacity = 0;
  }
}

static
void monomial_assign(const polynomial_context_t* ctx, monomial_t* m, const monomial_t* from, int sort) {
  if (m != from) {
    monomial_destruct(m);
    monomial_construct_copy(ctx, m, from, sort);
  }
}

static
void monomial_push(monomial_t* m, variable_t x, unsigned d) {
  if (m->n == m->capacity) {
    m->capacity += 5;
    m->p = realloc(m->p, m->capacity*sizeof(power_t));
  }
  assert(m->n < m->capacity);
  m->p[m->n].x = x;
  m->p[m->n].d = d;
  ++ (m->n);
}

static
void monomial_pop(monomial_t* m) {
  assert(m->n > 0);
  -- m->n;
}

static
int monomial_print(const polynomial_context_t* ctx, const monomial_t* m, FILE* out) {
  int ret = 0;
  ret += integer_print(&m->a, out);
  ret += fprintf(out, " * ");
  int i = 0;
  for (i = 0; i < m->n; ++ i) {
    if (i) {
      ret += fprintf(out, "*");
    }
    ret += fprintf(out, "%s^%u", variable_db_ops.get_name(ctx->var_db, m->p[i].x), m->p[i].d);
  }
  return ret;
}

#define SWAP(m1, m2) { monomial_t tmp = m1; m1 = m2; m2 = tmp; }
#define MIN(x, y) (x < y ? x : y)

void monomial_gcd(const polynomial_context_t* ctx, monomial_t* gcd, const monomial_t* m1, const monomial_t* m2) {

  assert(ctx->K == Z);

  monomial_t result;
  monomial_construct(ctx, &result);

  // GCD of the coefficients
  integer_gcd_Z(&result.a, &m1->a, &m2->a);

  // GCD of the power
  int m1_i = 0, m2_i = 0;
  while (m1_i < m1->n && m2_i < m2->n) {
    // Only keep powers that are equal
    // Variables in the monomial go top to bottom
    int var_cmp = ctx->var_order->ops->cmp(ctx->var_order, m1->p[m1_i].x, m2->p[m2_i].x);
    if (var_cmp == 0) {
      variable_t x = m1->p[m1_i].x;
      size_t d = MIN(m1->p[m1_i].d, m2->p[m2_i].d);
      monomial_push(&result, x, d);
      m1_i ++;
      m2_i ++;
    } else if (var_cmp > 0) {
      m1_i ++;
    } else if (var_cmp < 0) {
      m2_i ++;
    }
  }

  SWAP(result, *gcd);
  monomial_destruct(&result);
}

const monomial_ops_t monomial_ops = {
    monomial_construct,
    monomial_construct_copy,
    monomial_destruct,
    monomial_clear,
    monomial_assign,
    monomial_print,
    monomial_push,
    monomial_pop,
    monomial_gcd,
};

