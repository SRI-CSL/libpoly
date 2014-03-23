/*
 * polynomial_coefficient.c
 *
 *  Created on: Feb 4, 2014
 *      Author: dejan
 */

#include "polynomial/polynomial.h"
#include "polynomial/monomial.h"
#include "polynomial/coefficient.h"

#include "number/integer.h"
#include "interval/interval.h"

#include "utils/assignment.h"
#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <assert.h>
#include <malloc.h>
#include <string.h>

/**
 * Make sure that the coefficient has the given capacity for the given variable.
 * If the polynomial a constant or a polynomial in a smaller variable, we
 * transform it into a polynomial and adjust the coefficent.
 */
static void
coefficient_ensure_capacity(const polynomial_context_t* ctx, coefficient_t* C, variable_t x, size_t capacity);

/**
 * Assumes:
 * * Sub-coefficients have been normalized.
 *
 * Normalize the coefficient:
 * * If the coefficient is a polynomial with only the constant coefficient, it
 *   should be upgraded
 * * The highest coefficient should be non-zero
 */
static void
coefficient_normalize(const polynomial_context_t* ctx, coefficient_t* C);

int
coefficient_is_normalized(const polynomial_context_t* ctx, coefficient_t* C);

STAT_DECLARE(int, coefficient, construct);

static inline
void coefficient_construct(const polynomial_context_t* ctx, coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_construct()\n");
  STAT(coefficient, construct) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_from_int(ctx->K, &C->value.num, 0);
}

STAT_DECLARE(int, coefficient, construct_from_int);

static inline
void coefficient_construct_from_int(const polynomial_context_t* ctx, coefficient_t* C, long C_int) {
  TRACE("coefficient::internal", "coefficient_construct_from_int()\n");
  STAT(coefficient, construct_from_int) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_from_int(ctx->K, &C->value.num, C_int);
}

STAT_DECLARE(int, coefficient, construct_from_integer);

static inline
void coefficient_construct_from_integer(const polynomial_context_t* ctx, coefficient_t* C, const integer_t* C_integer) {
  TRACE("coefficient::internal", "coefficient_construct_from_integer()\n");
  STAT(coefficient, construct_from_integer) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_copy(ctx->K, &C->value.num, C_integer);
}

STAT_DECLARE(int, coefficient, construct_rec);

static
void coefficient_construct_rec(const polynomial_context_t* ctx, coefficient_t* C, variable_t x, size_t capacity) {
  TRACE("coefficient::internal", "coefficient_construct_rec()\n");
  STAT(coefficient, construct_rec) ++;

  C->type = COEFFICIENT_POLYNOMIAL;
  C->value.rec.x = x;
  C->value.rec.size = 0;
  C->value.rec.capacity = 0;
  C->value.rec.coefficients = 0;
  coefficient_ensure_capacity(ctx, C, x, capacity);
}

STAT_DECLARE(int, coefficient, construct_simple);

static
void coefficient_construct_simple(const polynomial_context_t* ctx, coefficient_t* C, const integer_t* a, variable_t x, unsigned n) {
  TRACE("coefficient::internal", "coefficient_construct_simple()\n");
  STAT(coefficient, construct_simple) ++;

  if (n == 0) {
    coefficient_construct_from_integer(ctx, C, a);
  } else {
    // x^n
    coefficient_construct_rec(ctx, C, x, n+1);
    integer_assign(ctx->K, &COEFF(C, n)->value.num, a);
  }
}

STAT_DECLARE(int, coefficient, construct_copy);

static
void coefficient_construct_copy(const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from) {
  TRACE("coefficient::internal", "coefficient_construct_copy()\n");
  STAT(coefficient, construct_copy) ++;

  int i;
  switch(from->type) {
  case COEFFICIENT_NUMERIC:
    C->type = COEFFICIENT_NUMERIC;
    integer_construct_copy(ctx->K, &C->value.num, &from->value.num);
    break;
  case COEFFICIENT_POLYNOMIAL:
    C->type = COEFFICIENT_POLYNOMIAL;
    C->value.rec.x = VAR(from);
    C->value.rec.size = SIZE(from);
    C->value.rec.capacity = SIZE(from);
    C->value.rec.coefficients = malloc(SIZE(from) * sizeof(coefficient_t));
    for (i = 0; i < SIZE(from); ++ i) {
      coefficient_construct_copy(ctx, COEFF(C, i), COEFF(from, i));
    }
    break;
  }
}

STAT_DECLARE(int, coefficient, construct_from_univariate);

static inline
void coefficient_construct_from_univariate(const polynomial_context_t* ctx,
    coefficient_t* C, const upolynomial_t* C_u, variable_t x) {

  TRACE("coefficient::internal", "coefficient_construct_from_int()\n");
  STAT(coefficient, construct_from_univariate) ++;

  // Get the coefficients
  size_t C_u_deg = upolynomial_ops.degree(C_u);
  integer_t* coeff = malloc(sizeof(integer_t)*(C_u_deg + 1));

  int i;
  for (i = 0; i <= C_u_deg; ++ i) {
    integer_construct_from_int(ctx->K, coeff + i, 0);
  }

  upolynomial_ops.unpack(C_u, coeff);

  // Construct the polynomial
  coefficient_construct_rec(ctx, C, x, C_u_deg + 1);

  // Move over the coefficients
  for (i = 0; i <= C_u_deg; ++ i) {
    integer_swap(ctx->K, &COEFF(C, i)->value.num, coeff + i);
    integer_destruct(coeff + i);
  }

  // Normalize (it might be constant)
  coefficient_normalize(ctx, C);

  assert(coefficient_is_normalized(ctx, C));
}

static
void coefficient_destruct(coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_destruct()\n");

  int i;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    integer_destruct(&C->value.num);
    break;
  case COEFFICIENT_POLYNOMIAL:
    for (i = 0; i < CAPACITY(C); ++ i) {
      coefficient_destruct(COEFF(C, i));
    }
    free(C->value.rec.coefficients);
  }
}

STAT_DECLARE(int, coefficient, swap);

static inline
void coefficient_swap(coefficient_t* C1, coefficient_t* C2) {
  TRACE("coefficient::internal", "coefficient_swap()\n");
  STAT(coefficient, swap) ++;
  coefficient_t tmp = *C1;
  *C1 = *C2;
  *C2 = tmp;
}

STAT_DECLARE(int, coefficient, assign);

static
void coefficient_assign(const polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from) {
  TRACE("coefficient::internal", "coefficient_assign()\n");
  STAT(coefficient, assign) ++;

  if (C != from) {
    coefficient_t result;
    switch(from->type) {
    case COEFFICIENT_NUMERIC:
      if (C->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(C);
        coefficient_construct_copy(ctx, C, from);
      } else {
        integer_assign(ctx->K, &C->value.num, &from->value.num);
      }
      break;
    case COEFFICIENT_POLYNOMIAL:
      coefficient_construct_copy(ctx, &result, from);
      coefficient_swap(&result, C);
      coefficient_destruct(&result);
      break;
    }
  }

  assert(coefficient_is_normalized(ctx, C));
}

STAT_DECLARE(int, coefficient, assign_int);

static
void coefficient_assign_int(const polynomial_context_t* ctx, coefficient_t* C, long x) {
  TRACE("coefficient::internal", "coefficient_assign_int()\n");
  STAT(coefficient, assign) ++;

  if (C->type == COEFFICIENT_POLYNOMIAL) {
    coefficient_destruct(C);
    coefficient_construct_from_int(ctx, C, x);
  } else {
    integer_assign_int(ctx->K, &C->value.num, x);
  }

  assert(coefficient_is_normalized(ctx, C));
}

const coefficient_t* coefficient_lc_safe(const polynomial_context_t* ctx, const coefficient_t* C, variable_t x) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return C;
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (VAR(C) == x) {
      return COEFF(C, SIZE(C) - 1);
    } else {
      assert(ctx->var_order->ops->cmp(ctx->var_order, x, VAR(C)) > 0);
      return C;
    }
  default:
    assert(0);
    return 0;
  }
}

const coefficient_t* coefficient_lc(const coefficient_t* C) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return C;
    break;
  case COEFFICIENT_POLYNOMIAL:
    return COEFF(C, SIZE(C) - 1);
    break;
  }
  assert(0);
  return 0;
}

static inline
int coefficient_is_constant(const coefficient_t* C) {
  return C->type == COEFFICIENT_NUMERIC;
}

static inline
size_t coefficient_degree(const coefficient_t* C) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 0;
    break;
  case COEFFICIENT_POLYNOMIAL:
    return SIZE(C) - 1;
    break;
  }
  assert(0);
  return 0;
}

static inline
size_t coefficient_degree_safe(const polynomial_context_t* ctx, const coefficient_t* C, variable_t x) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 0;
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (VAR(C) == x) {
      return SIZE(C) - 1;
    } else {
      assert(ctx->var_order->ops->cmp(ctx->var_order, x, VAR(C)) > 0);
      return 0;
    }
    break;
  }
  assert(0);
  return 0;
}

static inline
variable_t coefficient_top_variable(const coefficient_t* C) {
  assert(C->type == COEFFICIENT_POLYNOMIAL);
  return VAR(C);
}

static
const coefficient_t* coefficient_get_coefficient(const coefficient_t* C, size_t d) {

  assert(d <= coefficient_degree(C));

  switch(C->type) {
  case COEFFICIENT_NUMERIC:
    return C;
    break;
  case COEFFICIENT_POLYNOMIAL:
    return COEFF(C, d);
    break;
  }

  assert(0);
  return 0;
}

STAT_DECLARE(int, coefficient, is_zero);

static inline
int coefficient_is_zero(const polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, is_zero) ++;
  return C->type == COEFFICIENT_NUMERIC && integer_is_zero(ctx->K, &C->value.num);
}

STAT_DECLARE(int, coefficient, is_one);

static inline
int coefficient_is_one(const polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, is_one) ++;
  return C->type == COEFFICIENT_NUMERIC && integer_cmp_int(ctx->K, &C->value.num, 1) == 0;
}

void coefficient_value_approx(const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m, interval_t* value) {

  if (C->type == COEFFICIENT_NUMERIC) {
    interval_t result;
    interval_construct_from_integer(&result, &C->value.num, 0, &C->value.num, 0);
    interval_swap(value, &result);
    interval_destruct(&result);
  } else {
    interval_t result, tmp1, tmp2, x_value;
    interval_construct_zero(&result);
    interval_construct_zero(&tmp1);
    interval_construct_zero(&tmp2);
    interval_construct_zero(&x_value);
    // Get the value of x
    assignment_get_value_approx(m, VAR(C), &x_value);

    // We compute using powers, just an attemp to compute better. For example
    // if p = x^2 + x and x = [-1, 1] then
    //  a) x(x + 1) = [-1, 1]*[0, 2] = [-2, 2]
    //  b) x^2 + x = [0, 1] + [-1, 1] = [-1, 2]
    // we choose to do it b) way

    // Compute
    int i;
    for (i = 0; i < SIZE(C); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_value_approx(ctx, COEFF(C, i), m, &tmp1);
        interval_pow(&tmp2, &x_value, i);
        interval_mul(&tmp2, &tmp2, &tmp1);
        interval_add(&result, &result, &tmp2);
      }
    }

    interval_swap(&result, value);
    interval_destruct(&x_value);
    interval_destruct(&tmp1);
    interval_destruct(&tmp2);
    interval_destruct(&result);
  }
}

STAT_DECLARE(int, coefficient, sgn);

int coefficient_sgn(const polynomial_context_t* ctx, const coefficient_t* C, const assignment_t* m) {

  TRACE("coefficient::internal", "coefficient_sgn()\n");
  STAT(coefficient, sgn) ++;

  int sgn;

  if (C->type == COEFFICIENT_NUMERIC) {
    // For numberic coefficients we're done
    sgn = integer_sgn(ctx->K, &C->value.num);
  } else {
    assert(C->type == COEFFICIENT_POLYNOMIAL);

    // Approximate the value of C
    interval_t C_approx;
    interval_construct_zero(&C_approx);

    // Approximate the value
    coefficient_value_approx(ctx, C, m, &C_approx);

    // Safe to give the sign based on the interval bound
    sgn = interval_sgn(&C_approx);

    // Destruct temps
    interval_destruct(&C_approx);
  }


  return sgn;
}

STAT_DECLARE(int, coefficient, lc_sgn);

int coefficient_lc_sgn(const polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, lc_sgn) ++;

  while (C->type != COEFFICIENT_NUMERIC) {
    C = coefficient_lc(C);
  }

  return integer_sgn(ctx->K, &C->value.num);
}


STAT_DECLARE(int, coefficient, in_order);

static
int coefficient_in_order(const polynomial_context_t* ctx, const coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_in_order()\n");
  STAT(coefficient, in_order) ++;

  int i;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 1;
    break;
  case COEFFICIENT_POLYNOMIAL:
    // Check that the top is bigger than the top of coefficient and run
    // recursively
    for (i = 0; i < SIZE(C); ++ i) {
      const coefficient_t* C_i = COEFF(C, i);
      if (C_i->type == COEFFICIENT_POLYNOMIAL) {
        if (ctx->var_order->ops->cmp(ctx->var_order, VAR(C), VAR(C_i)) <= 0) {
          // Top variable must be bigger than others
          return 0;
        } else if (!coefficient_in_order(ctx, C_i)) {
          return 0;
        }
      }
    }

    break;
  }

  return 1;
}


static
int coefficient_cmp_general(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2, int compare_values) {
  int cmp;

  if (C1->type == COEFFICIENT_NUMERIC && C2->type == COEFFICIENT_NUMERIC) {
    if (compare_values) {
      cmp = integer_cmp(ctx->K, &C1->value.num, &C2->value.num);
    } else {
      cmp = 0;
    }
  } else if (C1->type == COEFFICIENT_NUMERIC) {
    // C1 is a constant, C1 is always smaller
    return -1;
  } else if (C2->type == COEFFICIENT_NUMERIC) {
    // C2 is a constant, C1 is always bigger
    return 1;
  } else {
    // Both are polynomials, compare the variable
    int var_cmp = ctx->var_order->ops->cmp(ctx->var_order, VAR(C1), VAR(C2));
    if (var_cmp == 0) {
      if (compare_values) {
        // If the variables are the same, compare lexicographically
        int deg_cmp = ((int) SIZE(C1)) - ((int) SIZE(C2));
        if (deg_cmp == 0) {
          int i = C1->value.rec.size - 1;
          for (; i >= 0; -- i) {
            int coeff_cmp = coefficient_cmp_general(ctx, COEFF(C1, i), COEFF(C2, i), compare_values);
            if (coeff_cmp != 0) {
              cmp = coeff_cmp;
              break;
            }
          }
          if (i < 0) {
            // All coefficients equal
            cmp = 0;
          }
        } else {
          cmp = deg_cmp;
        }
      } else {
        return 0;
      }
    } else {
      // Variable comparison is enough
      cmp = var_cmp;
    }
  }

  TRACE("coefficien::internal", "coefficient_cmp() => %d\n", cmp);
  return cmp;
}

STAT_DECLARE(int, coefficient, cmp);

static inline
int coefficient_cmp(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_cmp()\n");
  STAT(coefficient, cmp) ++;

  return coefficient_cmp_general(ctx, C1, C2, 1);
}

STAT_DECLARE(int, coefficient, cmp_type);

static inline
int coefficient_cmp_type(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::internal", "coefficient_cmp_type()\n");
  STAT(coefficient, cmp_type) ++;

  return coefficient_cmp_general(ctx, C1, C2, 0);
}

STAT_DECLARE(int, coefficient, divides);

static
int coefficient_divides(const polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_divides()\n");
  STAT(coefficient, divides) ++;

  int divides = 0;

  return divides;
}

static
char* power_symbol = 0;

/** Set the power symbol for print-outs */
static
void coefficient_set_power_symbol(const char* pow) {
  if (power_symbol) {
    free(power_symbol);
    power_symbol = 0;
  }
  if (pow) {
    power_symbol = strdup(pow);
  }
}

static
const char* coefficient_get_power_symbol(void) {
  if (power_symbol) return power_symbol;
  else return "^";
}

static
int coefficient_print(const polynomial_context_t* ctx, const coefficient_t* C, FILE* out) {
  int i, k = 0, ret = 0;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    ret += integer_print(&C->value.num, out);
    break;
  case COEFFICIENT_POLYNOMIAL: {
    // The polynomial
    const char* var_name = variable_db_ops.get_name(ctx->var_db, C->value.rec.x);
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
              integer_t tmp;
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
            ret += fprintf(out, "*%s%s%d", var_name, coefficient_get_power_symbol(), i);
          }
        }
      }
    }
    break;
  }
  }
  return ret;
}

static
char* coefficient_to_string(const polynomial_context_t* ctx, const coefficient_t* C) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  coefficient_print(ctx, C, f);
  fclose(f);
  return str;
}

/** Function type called on coefficient traversal */
typedef void (*traverse_f) (const polynomial_context_t* ctx, monomial_t* p, void* data);

static
void coefficient_traverse(const polynomial_context_t* ctx, const coefficient_t* C, traverse_f f, monomial_t* m, void* data) {

  if (debug_trace_ops.is_enabled("coefficient::order")) {
    tracef("order = "); variable_order_simple_ops.print((variable_order_simple_t*) ctx->var_order, ctx->var_db, trace_out); tracef("\n");
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
    tracef("m = "); monomial_ops.print(ctx, m, trace_out); tracef("\n");
  }

  int d;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    integer_assign(ctx->K, &m->a, &C->value.num);
    (*f)(ctx, m, data);
    break;
  case COEFFICIENT_POLYNOMIAL:
    // The constant
    if (!coefficient_is_zero(ctx, COEFF(C, 0))) {
      coefficient_traverse(ctx, COEFF(C, 0), f, m, data);
    }
    // Power of x
    for (d = 1; d < SIZE(C); ++ d) {
      if (!coefficient_is_zero(ctx, COEFF(C, d))) {
        monomial_ops.push(m, VAR(C), d);
        coefficient_traverse(ctx, COEFF(C, d), f, m, data);
        monomial_ops.pop(m);
      }
    }
    break;
  }
}

/**
 * Method called to add a monomial to C. The monomial should be ordered in the
 * same order as C, top variable at the m[0].
 */
static
void coefficient_add_monomial(const polynomial_context_t* ctx, monomial_t* m, void* C_void) {

  coefficient_t* C = (coefficient_t*) C_void;

  if (debug_trace_ops.is_enabled("coefficient::order")) {
    tracef("coefficient_add_monomial():\n");
    tracef("m = "); monomial_ops.print(ctx, m, trace_out); tracef("\n");
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
  }

  if (m->n == 0) {
    // Just add a constant to C
    switch (C->type) {
    case COEFFICIENT_NUMERIC:
      integer_add(ctx->K, &C->value.num, &C->value.num, &m->a);
      break;
    case COEFFICIENT_POLYNOMIAL:
      coefficient_add_monomial(ctx, m, COEFF(C, 0));
      break;
    }
  } else {
    // Proper monomial (first variable is the top one)
    variable_t x = m->p[0].x;
    unsigned d = m->p[0].d;
    // Compare the variables
    if (C->type == COEFFICIENT_NUMERIC || ctx->var_order->ops->cmp(ctx->var_order, x, VAR(C)) >= 0) {
      coefficient_ensure_capacity(ctx, C, x, d+1);
      // Now, add the monomial to the right place
      m->p ++;
      m->n --;
      coefficient_add_monomial(ctx, m, COEFF(C, d));
      coefficient_normalize(ctx, C);
      m->p --;
      m->n ++;
    } else {
      coefficient_add_monomial(ctx, m, COEFF(C, 0));
    }
  }

  assert(coefficient_is_normalized(ctx, C));
}

static
void coefficient_order_and_add_monomial(const polynomial_context_t* ctx, monomial_t* m, void* C_void) {
  monomial_t m_ordered;
  monomial_ops.construct_copy(ctx, &m_ordered, m, /** sort */ 1);
  coefficient_add_monomial(ctx, &m_ordered, C_void);
  monomial_ops.destruct(&m_ordered);
}

STAT_DECLARE(int, coefficient, order);

static
void coefficient_order(const polynomial_context_t* ctx, coefficient_t* C) {
  TRACE("coefficient", "coefficient_order()\n");
  STAT(coefficient, order) ++;

  if (C->type == COEFFICIENT_NUMERIC) {
    // Numeric coefficients are always OK
    return;
  }

  if (debug_trace_ops.is_enabled("coefficient::order")) {
    tracef("order = "); variable_order_simple_ops.print((variable_order_simple_t*) ctx->var_order, ctx->var_db, trace_out); tracef("\n");
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
  }

  // The coefficient we are building
  coefficient_t result;
  coefficient_construct(ctx, &result);
  // The monomials build in the original order
  monomial_t m_tmp;
  monomial_ops.construct(ctx, &m_tmp);
  // For each monomial of C, add it to the result
  coefficient_traverse(ctx, C, coefficient_order_and_add_monomial, &m_tmp, &result);
  // Keep the result
  coefficient_swap(C, &result);
  // Destroy temps
  monomial_ops.destruct(&m_tmp);
  coefficient_destruct(&result);

  assert(coefficient_is_normalized(ctx, C));
}

STAT_DECLARE(int, coefficient, add);

#define MAX(x, y) (x >= y ? x : y)

static
void coefficient_add(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_add()\n");
  STAT(coefficient, add) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  coefficient_t result;

  int type_cmp = coefficient_cmp_type(ctx, C1, C2);

  if (type_cmp == 0) {
    if (C1->type == COEFFICIENT_NUMERIC) {
      assert(C2->type == COEFFICIENT_NUMERIC);
      // Add the integers
      integer_add(ctx->K, &S->value.num, &C1->value.num, &C2->value.num);
    } else {
      assert(C1->type == COEFFICIENT_POLYNOMIAL);
      assert(C2->type == COEFFICIENT_POLYNOMIAL);
      assert(VAR(C1) == VAR(C2));
      // Two polynomials over the same top variable
      size_t max_size = MAX(SIZE(C1), SIZE(C2));
      coefficient_construct_rec(ctx, &result, VAR(C1), max_size);
      int i;
      for (i = 0; i < max_size; ++ i) {
        if (i < SIZE(C1)) {
          if (i < SIZE(C2)) {
            // add C1 and C2
            coefficient_add(ctx, COEFF(&result, i), COEFF(C1, i), COEFF(C2, i));
          } else {
            // copy C1 coefficient
            coefficient_assign(ctx, COEFF(&result, i), COEFF(C1, i));
          }
        } else {
          // copy C2 coefficient
          coefficient_assign(ctx, COEFF(&result, i), COEFF(C2, i));
        }
      }
      coefficient_normalize(ctx, &result);
      coefficient_swap(&result, S);
      coefficient_destruct(&result);
    }
  } else if (type_cmp > 0) {
    // C1 > C2, add C2 into the constant of C1
    // We can't assign S to C1, since C2 might be S, so we use a temp
    coefficient_construct_copy(ctx, &result, C1);
    coefficient_add(ctx, COEFF(&result, 0), COEFF(C1, 0), C2);
    coefficient_swap(&result, S);
    coefficient_destruct(&result);
    // Since C1 is not a constant, no normalization needed, same size
  } else {
    // C1 < C2, add C1 into the constant of C2
    // We can't assign C2 to S1, since C1 might be S, so we use a temp
    coefficient_construct_copy(ctx, &result, C2);
    coefficient_add(ctx, COEFF(&result, 0), C1, COEFF(C2, 0));
    coefficient_swap(&result, S);
    coefficient_destruct(&result);
    // Since C2 is not a constant, no normalization needed, same size
  }

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("add = "); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, neg);

static
void coefficient_neg(const polynomial_context_t* ctx, coefficient_t* N, const coefficient_t* C) {
  TRACE("coefficient::arith", "coefficient_neg()\n");
  STAT(coefficient, neg) ++;

  int i;
  coefficient_t result;

  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    if (N->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(N);
      coefficient_construct(ctx, N);
    }
    integer_neg(ctx->K, &N->value.num, &C->value.num);
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (N != C) {
      coefficient_construct_rec(ctx, &result, VAR(C), SIZE(C));
      for (i = 0; i < SIZE(C); ++i) {
        if (!coefficient_is_zero(ctx, COEFF(C, i))) {
          coefficient_neg(ctx, COEFF(&result, i), COEFF(C, i));
        }
      }
      coefficient_normalize(ctx, &result);
      coefficient_swap(&result, N);
      coefficient_destruct(&result);
    } else {
      // In-place negation
      for (i = 0; i < SIZE(C); ++i) {
        if (!coefficient_is_zero(ctx, COEFF(C, i))) {
          coefficient_neg(ctx, COEFF(N, i), COEFF(C, i));
        }
      }
    }
    break;
  }

  assert(coefficient_is_normalized(ctx, N));
}

STAT_DECLARE(int, coefficient, sub);

static
void coefficient_sub(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_sub()\n");
  STAT(coefficient, sub) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  coefficient_t result;

  int type_cmp = coefficient_cmp_type(ctx, C1, C2);

  if (type_cmp == 0) {
    if (C1->type == COEFFICIENT_NUMERIC) {
      assert(C2->type == COEFFICIENT_NUMERIC);
      // Subtract the integers
      integer_sub(ctx->K, &S->value.num, &C1->value.num, &C2->value.num);
    } else {
      assert(C1->type == COEFFICIENT_POLYNOMIAL);
      assert(C2->type == COEFFICIENT_POLYNOMIAL);
      // Two polynomials over the same top variable
      assert(VAR(C1) == VAR(C2));
      size_t max_size = MAX(SIZE(C1), SIZE(C2));
      coefficient_construct_rec(ctx, &result, VAR(C1), max_size);
      int i;
      for (i = 0; i < max_size; ++ i) {
        if (i < SIZE(C1)) {
          if (i < SIZE(C2)) {
            // subtract C1 and C2
            coefficient_sub(ctx, COEFF(&result, i), COEFF(C1, i), COEFF(C2, i));
          } else {
            // copy C1 coefficient
            coefficient_assign(ctx, COEFF(&result, i), COEFF(C1, i));
          }
        } else {
          // copy -C2 ceofficient
          coefficient_neg(ctx, COEFF(&result, i), COEFF(C2, i));
        }
      }
      coefficient_normalize(ctx, &result);
      coefficient_swap(&result, S);
      coefficient_destruct(&result);
    }
  } else if (type_cmp > 0) {
    // C1 > C2, subtract C2 into the constant of C1
    // Can't assign C1 to S, since C2 might be S
    coefficient_construct_copy(ctx, &result, C1);
    coefficient_sub(ctx, COEFF(&result, 0), COEFF(C1, 0), C2);
    coefficient_swap(&result, S);
    coefficient_destruct(&result);
    // Since C1 is not a constant, no normalization is needed
  } else {
    // C1 < C2, compute -(C2 - C1)
    coefficient_sub(ctx, S, C2, C1);
    coefficient_neg(ctx, S, S);
    // Since C2 is not a constant, no normalization is needed
  }

  assert(coefficient_is_normalized(ctx, S));
}


static
void coefficient_add_mul(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

STAT_DECLARE(int, coefficient, mul);

static
void coefficient_mul(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_mul()\n");
  STAT(coefficient, mul) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  int i, j;
  coefficient_t result;

  int type_cmp = coefficient_cmp_type(ctx, C1, C2);

  if (type_cmp == 0) {
    if (C1->type == COEFFICIENT_NUMERIC) {
      assert(C2->type == COEFFICIENT_NUMERIC);
      // Multiply the integers
      integer_mul(ctx->K, &P->value.num, &C1->value.num, &C2->value.num);
    } else {
      assert(C1->type == COEFFICIENT_POLYNOMIAL);
      assert(C2->type == COEFFICIENT_POLYNOMIAL);
      // Two polynomials over the same top variable
      assert(VAR(C1) == VAR(C2));
      coefficient_construct_rec(ctx, &result, VAR(C1), SIZE(C1) + SIZE(C2) - 1);
      for (i = 0; i < SIZE(C1); ++ i) {
        if (!coefficient_is_zero(ctx, COEFF(C1, i))) {
          for (j = 0; j < SIZE(C2); ++ j) {
            if (!coefficient_is_zero(ctx, COEFF(C2, j))) {
              coefficient_add_mul(ctx, COEFF(&result, i + j), COEFF(C1, i), COEFF(C2, j));
              if (debug_trace_ops.is_enabled("coefficient::arith")) {
                tracef("result = "); coefficient_print(ctx, &result, trace_out); tracef("\n");
              }
            }
          }
        }
      }
      coefficient_normalize(ctx, &result);
      coefficient_swap(&result, P);
      coefficient_destruct(&result);
    }
  } else if (type_cmp > 0) {
    assert(C1->type == COEFFICIENT_POLYNOMIAL);
    // C1 > C2, multiply each coefficient of C1 with C2
    coefficient_construct_rec(ctx, &result, VAR(C1), SIZE(C1));
    for (i = 0; i < SIZE(C1); ++ i) {
      coefficient_mul(ctx, COEFF(&result, i), COEFF(C1, i), C2);
    }
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, P);
    coefficient_destruct(&result);
  } else {
    // C1 < C2, multuply each coefficient of C2 with C1
    coefficient_construct_rec(ctx, &result, VAR(C2), SIZE(C2));
    for (i = 0; i < SIZE(C2); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C2, i))) {
        coefficient_mul(ctx, COEFF(&result, i), C1, COEFF(C2, i));
      }
    }
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, P);
    coefficient_destruct(&result);
  }

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("mul = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, P));
}

STAT_DECLARE(int, coefficient, mul_int);

static
void coefficient_mul_int(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, long a) {
  TRACE("coefficient::arith", "coefficient_mul_int()\n");
  STAT(coefficient, mul_int) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
    tracef("n  = %ld", a);
  }

  int i;
  coefficient_t result;

  if (C->type == COEFFICIENT_NUMERIC) {
    if (P->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_construct(ctx, &result);
      integer_mul_int(ctx->K, &result.value.num, &C->value.num, a);
      coefficient_swap(&result, P);
      coefficient_destruct(&result);
    } else {
      integer_mul_int(ctx->K, &P->value.num, &C->value.num, a);
    }
  } else {
    coefficient_construct_rec(ctx, &result, VAR(C), SIZE(C));
    for (i = 0; i < SIZE(C); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_mul_int(ctx, COEFF(&result, i), COEFF(C, i), a);
      }
    }
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, P);
    coefficient_destruct(&result);
  }

  assert(coefficient_is_normalized(ctx, P));
}


STAT_DECLARE(int, coefficient, shl);

static
void coefficient_shl(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, variable_t x, unsigned n) {
  TRACE("coefficient::arith", "coefficient_shl()\n");
  STAT(coefficient, shl) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
    tracef("x = %s\n", variable_db_ops.get_name(ctx->var_db, x));
    tracef("n  = %u\n", n);
  }

  coefficient_assign(ctx, S, C);
  if (n > 0) {
    int old_size = (S->type == COEFFICIENT_NUMERIC || VAR(S) != x) ? 1 : SIZE(S);
    coefficient_ensure_capacity(ctx, S, x, old_size + n);
    int i;
    for (i = old_size - 1; i >= 0; -- i) {
      if (!coefficient_is_zero(ctx, COEFF(S, i))) {
        coefficient_swap(COEFF(S, i + n), COEFF(S, i));
      }
    }
  }

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("coefficient_shl() =>"); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, shr);

static
void coefficient_shr(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, unsigned n) {
  TRACE("coefficient::arith", "coefficient_shr()\n");
  STAT(coefficient, shl) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
    tracef("n  = %u\n", n);
  }

  assert(C->type == COEFFICIENT_POLYNOMIAL);
  assert(n < SIZE(C));

  if (n == 0) {
    coefficient_assign(ctx, S, C);
  } else if (n + 1 == SIZE(C)) {
    if (S == C) {
      coefficient_t result;
      coefficient_construct_copy(ctx, &result, coefficient_lc(C));
      coefficient_swap(&result, S);
      coefficient_destruct(&result);
    } else {
      coefficient_assign(ctx, S, coefficient_lc(C));
    }
  } else {
    coefficient_t result;
    coefficient_construct_rec(ctx, &result, VAR(C), SIZE(C) - n);
    int i;
    for (i = 0; i < SIZE(C) - n; ++ i) {
      coefficient_assign(ctx, COEFF(&result, i), COEFF(C, i + n));
    }
    coefficient_swap(&result, S);
    coefficient_destruct(&result);
  }

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("coefficient_shr() =>"); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, pow);

static
void coefficient_pow(const polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, unsigned n) {
  TRACE("coefficient", "coefficient_pow()\n");
  STAT(coefficient, pow) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
  }

  if (n == 0) {
    coefficient_assign_int(ctx, P, 1);
    return;
  }

  if (n == 1) {
    coefficient_assign(ctx, P, C);
    return;
  }

  coefficient_t result, tmp;

  switch(C->type) {
  case COEFFICIENT_NUMERIC:
    if (P->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_construct(ctx, &result);
      integer_pow(ctx->K, &result.value.num, &C->value.num, n);
      coefficient_swap(P, &result);
      coefficient_destruct(&result);
    } else {
      integer_pow(ctx->K, &P->value.num, &C->value.num, n);
    }
    break;
  case COEFFICIENT_POLYNOMIAL:
    // Accumulator for C^n (start with 1)
    coefficient_construct_from_int(ctx, &result, 1);
    coefficient_ensure_capacity(ctx, &result, VAR(C), (SIZE(C)-1)*n + 1);
    // C^power of 2 (start with C)
    coefficient_construct_copy(ctx, &tmp, C);
    while (n) {
      if (n & 1) {
        coefficient_mul(ctx, &result, &result, &tmp);
      }
      coefficient_mul(ctx, &tmp, &tmp, &tmp);
      n >>= 1;
    }
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, P);
    coefficient_destruct(&tmp);
    coefficient_destruct(&result);
    break;
  }

  assert(coefficient_is_normalized(ctx, P));
}

STAT_DECLARE(int, coefficient, add_mul);

static
void coefficient_add_mul(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_add_mul()\n");
  STAT(coefficient, add_mul) ++;

  if (debug_trace_ops.is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_ops.print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  if (C1->type == COEFFICIENT_NUMERIC && C2->type == COEFFICIENT_NUMERIC && S->type == COEFFICIENT_NUMERIC) {
    integer_add_mul(ctx->K, &S->value.num, &C1->value.num, &C2->value.num);
  } else {
    coefficient_t mul;
    coefficient_construct(ctx, &mul);
    coefficient_mul(ctx, &mul, C1, C2);
    coefficient_add(ctx, S, S, &mul);
    coefficient_destruct(&mul);
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, sub_mul);

static
void coefficient_sub_mul(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_sub_mul()\n");
  STAT(coefficient, sub_mul) ++;

  if (S->type == COEFFICIENT_NUMERIC && C1->type == COEFFICIENT_NUMERIC && C2->type == COEFFICIENT_NUMERIC) {
    integer_sub_mul(ctx->K, &S->value.num, &C1->value.num, &C2->value.num);
  } else {
    coefficient_t mul;
    coefficient_construct(ctx, &mul);
    coefficient_mul(ctx, &mul, C1, C2);
    coefficient_sub(ctx, S, S, &mul);
    coefficient_destruct(&mul);
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, derivative);

void coefficient_derivative(const polynomial_context_t* ctx, coefficient_t* C_d, const coefficient_t* C) {
  TRACE("coefficient", "coefficient_derivative()\n");
  STAT(coefficient, derivative) ++;

  int i;
  coefficient_t result;

  switch(C->type) {
  case COEFFICIENT_NUMERIC:
    coefficient_construct(ctx, &result);
    break;
  case COEFFICIENT_POLYNOMIAL:
    coefficient_construct_rec(ctx, &result, VAR(C), SIZE(C));
    for (i = 1; i < SIZE(C); ++ i) {
      coefficient_mul_int(ctx, COEFF(&result, i-1), COEFF(C, i), i);
    }
    coefficient_normalize(ctx, &result);
    break;
  }

  coefficient_swap(C_d, &result);
  coefficient_destruct(&result);

  assert(coefficient_is_normalized(ctx, C_d));
}

///
/// Forward declarations of division/reduction/gcd stuff
///

static
void coefficient_div(const polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2);

static
void coefficient_pp_cont(const polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C);

static
void coefficient_pp(const polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C);

static
void coefficient_cont(const polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C);

static
void coefficient_gcd(const polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2);

static
void coefficient_lcm(const polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2);

//
// Implementation of the division/reduction/gcd stuff
//

STAT_DECLARE(int, coefficient, reduce);

static
void coefficient_reduce(
    const polynomial_context_t* ctx,
    const coefficient_t* A, const coefficient_t* B,
    coefficient_t* P, coefficient_t* Q, coefficient_t* R,
    remaindering_type_t type)
{
  TRACE("coefficient", "coefficient_reduce()\n");
  STAT(coefficient, reduce) ++;

  if (debug_trace_ops.is_enabled("coefficient::reduce")) {
    tracef("A = "); coefficient_ops.print(ctx, A, trace_out); tracef("\n");
    tracef("B = "); coefficient_ops.print(ctx, B, trace_out); tracef("\n");
  }

  assert(A->type == COEFFICIENT_POLYNOMIAL);

  // Start with R = A, P = 1
  coefficient_t P_tmp, Q_tmp, R_tmp;
  if (P) {
    coefficient_construct_from_int(ctx, &P_tmp, 1);
  }
  if (Q) {
    coefficient_construct(ctx, &Q_tmp);
  }
  coefficient_construct_copy(ctx, &R_tmp, A);

  // The main variable we are reducing
  variable_t x = VAR(A);

  // Keep track of the degrees of the reduct to account for dense/sparse
  int R_deg = coefficient_degree(&R_tmp);
  int R_deg_prev = R_deg;
  int B_deg = coefficient_degree_safe(ctx, B, x);

  // Temporaries for computation
  coefficient_t lcm, r, b;
  coefficient_construct(ctx, &lcm);
  coefficient_construct_from_int(ctx, &r, 1);
  coefficient_construct(ctx, &b);

  // Invariant:
  //
  // P*A = Q*B + R
  //
  // initially
  //
  // P = 1, Q = 0, R = A
  //
  do {

    // Leading coefficient of B
    const coefficient_t* lc_B = coefficient_lc_safe(ctx, B, x);

    // Account for the sparse operation
    switch (type) {
    case REMAINDERING_PSEUDO_DENSE:
      if (R_deg_prev - R_deg > 1) {
        // Multiply with the missed power of lc(B)
        int missed = R_deg < B_deg ?
            R_deg_prev - B_deg : R_deg_prev - R_deg - 1;
        assert(missed >= 0);
        if (missed > 0) {
          coefficient_t pow;
          coefficient_construct(ctx, &pow);
          coefficient_pow(ctx, &pow, lc_B, missed);
          if (P) {
            coefficient_mul(ctx, &P_tmp, &P_tmp, &pow);
          }
          if (Q) {
            coefficient_mul(ctx, &Q_tmp, &Q_tmp, &pow);
          }
          coefficient_mul(ctx, &R_tmp, &R_tmp, &pow);
          coefficient_destruct(&pow);
        }
      }
      break;
    default:
      break;
    }

    // If we eliminated all of x we are done
    if (coefficient_is_zero(ctx, &R_tmp) || R_deg < B_deg) {
      break;
    }

    // How much x do we need to supply
    int d = R_deg - B_deg;

    // Eliminate the coefficient. We have
    //
    //   R = r_i * x^i + red(R)
    //   B = b_j * x^j + ref(B)
    //
    // and we compute
    //
    //   lcm = lcm(r_i, b_j)
    //   r = lcm/r_i
    //   b = (lcm/b_j)*x^d
    //
    //   R' = r*R - b*x^d*B
    //
    // thus eliminating the i-th power. As a consequence, the invariant is
    //
    //   P*A = Q*B + R        [*b]
    //
    //   r*P*A = r*Q*B + rR = r*Q*B + R' + b*B
    //   r*P*A = (r*Q + b)B + R'
    //
    // so we set (making sure that b is positive)
    //
    //   P' = b*P
    //   Q' = b*Q + c
    //

    // Leading coefficient of R
    const coefficient_t* lc_R = coefficient_lc_safe(ctx, &R_tmp, x);

    switch (type) {
    case REMAINDERING_EXACT_SPARSE:
    {
      // If we are exact, we assume that b_j divides r_i
      coefficient_div(ctx, &b, lc_R, lc_B);
      break;
    }
    case REMAINDERING_LCM_SPARSE:
    {
      // a, b, c
      coefficient_lcm(ctx, &lcm, lc_R, lc_B);
      coefficient_div(ctx, &r, &lcm, lc_R);
      coefficient_div(ctx, &b, &lcm, lc_B);
      if (coefficient_lc_sgn(ctx, &r) < 0) {
        coefficient_neg(ctx, &r, &r);
        coefficient_neg(ctx, &b, &b);
      }
      break;
    }
    case REMAINDERING_PSEUDO_DENSE:
    case REMAINDERING_PSEUDO_SPARSE:
    {
      coefficient_assign(ctx, &r, lc_B);
      coefficient_assign(ctx, &b, lc_R);
      break;
    }
    default:
      assert(0);
    }

    // Add the power of x
    coefficient_shl(ctx, &b, &b, x, d);

    assert(!coefficient_is_zero(ctx, &r));
    assert(!coefficient_is_zero(ctx, &b));

    if (debug_trace_ops.is_enabled("coefficient::reduce")) {
      tracef("lcm = "); coefficient_print(ctx, &lcm, trace_out); tracef("\n");
      tracef("R = "); coefficient_print(ctx, &R_tmp, trace_out); tracef("\n");
      tracef("r = "); coefficient_print(ctx, &r, trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, B, trace_out); tracef("\n");
      tracef("b = "); coefficient_print(ctx, &b, trace_out); tracef("\n");
    }

    // R'
    coefficient_mul(ctx, &R_tmp, &R_tmp, &r);
    coefficient_sub_mul(ctx, &R_tmp, B, &b);

    if (debug_trace_ops.is_enabled("coefficient::reduce")) {
      tracef("R' = "); coefficient_print(ctx, &R_tmp, trace_out); tracef("\n");
    }

    // Update the degrees of R
    R_deg_prev = R_deg;
    R_deg = coefficient_degree_safe(ctx, &R_tmp, x);
    assert(coefficient_is_zero(ctx, &R_tmp) || R_deg < R_deg_prev);


    // P' (if needed)
    if (P) {
      coefficient_mul(ctx, &P_tmp, &P_tmp, &r);
    }

    // Q' (if needed)
    if (Q) {
      coefficient_mul(ctx, &Q_tmp, &Q_tmp, &r);
      coefficient_add(ctx, &Q_tmp, &Q_tmp, &b);
    }

    if (debug_trace_ops.is_enabled("coefficient::reduce")) {
      if (P) {
        tracef("P = "); coefficient_print(ctx, &P_tmp, trace_out); tracef("\n");
      }
      if (Q) {
        tracef("Q = "); coefficient_print(ctx, &Q_tmp, trace_out); tracef("\n");
      }
    }

  } while (1);

  // Move the result out, and remove temps
  if (P) {
    coefficient_swap(P, &P_tmp);
    coefficient_destruct(&P_tmp);
    assert(coefficient_is_normalized(ctx, P));
  }
  if (Q) {
    coefficient_swap(Q, &Q_tmp);
    coefficient_destruct(&Q_tmp);
    assert(coefficient_is_normalized(ctx, Q));
  }
  if (R) {
    coefficient_swap(R, &R_tmp);
    assert(coefficient_is_normalized(ctx, R));
  }
  coefficient_destruct(&R_tmp);

  coefficient_destruct(&lcm);
  coefficient_destruct(&r);
  coefficient_destruct(&b);
}


STAT_DECLARE(int, coefficient, div);

static
void coefficient_div(const polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_div()\n");
  STAT(coefficient, div) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(!coefficient_is_zero(ctx, C2));

  if (coefficient_is_one(ctx, C2) || coefficient_is_zero(ctx, C1)) {
    coefficient_assign(ctx, D, C1);
    return;
  }

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  assert(cmp_type >= 0);

  if (cmp_type == 0 && C1->type == COEFFICIENT_NUMERIC) {
    assert(C2->type == COEFFICIENT_NUMERIC);
    if (D->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(D);
      coefficient_construct(ctx, D);
    }
    assert(integer_divides(ctx->K, &C2->value.num, &C1->value.num));
    integer_div_exact(ctx->K, &D->value.num, &C1->value.num, &C2->value.num);
  } else {
    if (debug_trace_ops.is_enabled("coefficient::check_division")) {
      coefficient_t R;
      coefficient_construct(ctx, &R);
      coefficient_reduce(ctx, C1, C2, 0, D, &R, REMAINDERING_EXACT_SPARSE);
      if (!coefficient_is_zero(ctx, &R)) {
        tracef("WRONG DIVISION!\n");
        tracef("P = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
        tracef("Q = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
      }
      coefficient_destruct(&R);
    } else {
      coefficient_reduce(ctx, C1, C2, 0, D, 0, REMAINDERING_EXACT_SPARSE);
    }
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_div() => "); coefficient_ops.print(ctx, D, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, D));
}

STAT_DECLARE(int, coefficient, rem);

static
void coefficient_rem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_rem()\n");
  STAT(coefficient, rem) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(!coefficient_is_zero(ctx, C2));

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  assert(cmp_type >= 0);

  if (cmp_type == 0 && C1->type == COEFFICIENT_NUMERIC) {
    assert(C2->type == COEFFICIENT_NUMERIC);
    if (R->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(R);
      coefficient_construct(ctx, R);
    }
    integer_rem_Z(&R->value.num, &C1->value.num, &C2->value.num);
  } else {
    coefficient_reduce(ctx, C1, C2, 0, 0, R, REMAINDERING_EXACT_SPARSE);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_rem() => "); coefficient_ops.print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, sprem);

static
void coefficient_sprem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_sprem()\n");
  STAT(coefficient, sprem) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(!coefficient_is_zero(ctx, C2));

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  assert(cmp_type >= 0);

  if (cmp_type == 0 && C1->type == COEFFICIENT_NUMERIC) {
    assert(C2->type == COEFFICIENT_NUMERIC);
    if (R->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(R);
      coefficient_construct(ctx, R);
    }
    integer_rem_Z(&R->value.num, &C1->value.num, &C2->value.num);
  } else {
    coefficient_reduce(ctx, C1, C2, 0, 0, R, REMAINDERING_PSEUDO_SPARSE);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_sprem() => "); coefficient_ops.print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, prem);

static
void coefficient_prem(const polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_prem()\n");
  STAT(coefficient, rem) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(!coefficient_is_zero(ctx, C2));

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  assert(cmp_type >= 0);

  if (cmp_type == 0 && C1->type == COEFFICIENT_NUMERIC) {
    assert(C2->type == COEFFICIENT_NUMERIC);
    if (R->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(R);
      coefficient_construct(ctx, R);
    }
    integer_rem_Z(&R->value.num, &C1->value.num, &C2->value.num);
  } else {
    coefficient_reduce(ctx, C1, C2, 0, 0, R, REMAINDERING_PSEUDO_DENSE);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_prem() => "); coefficient_ops.print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}


STAT_DECLARE(int, coefficient, divrem);

static
void coefficient_divrem(const polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_divrem()\n");
  STAT(coefficient, divrem) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(!coefficient_is_zero(ctx, C2));

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  assert(cmp_type >= 0);

  if (cmp_type == 0) {
    switch(C1->type) {
    case COEFFICIENT_NUMERIC:
      assert(C2->type == COEFFICIENT_NUMERIC);
      if (R->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(R);
        coefficient_construct(ctx, R);
      }
      integer_rem_Z(&R->value.num, &C1->value.num, &C2->value.num);
      break;
    case COEFFICIENT_POLYNOMIAL:
    {
      coefficient_reduce(ctx, C1, C2, 0, D, R, REMAINDERING_EXACT_SPARSE);
      break;
    }
    default:
      assert(0);
    }
  } else {
    // Just use the regular methods
    coefficient_rem(ctx, R, COEFF(C1, 0), C2);
    coefficient_div(ctx, D, C1, C2);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_divrem() => \n");
    tracef("D = "); coefficient_ops.print(ctx, D, trace_out); tracef("\n");
    tracef("D = "); coefficient_ops.print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, pp_cont);

static
void coefficient_pp_cont(const polynomial_context_t* ctx, coefficient_t* pp, coefficient_t* cont, const coefficient_t* C) {

  TRACE("coefficient", "coefficient_pp_cont()\n");
  STAT(coefficient, pp_cont) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C = "); coefficient_ops.print(ctx, C, trace_out); tracef("\n");
  }

  assert(ctx->K == Z);

  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    if (cont) {
      if (cont->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(cont);
        coefficient_construct_copy(ctx, cont, C);
      } else {
        integer_assign(ctx->K, &cont->value.num, &C->value.num);
      }
    }
    if (pp) {
      if (pp->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(pp);
        coefficient_construct_from_int(ctx, pp, 1);
      } else {
        integer_assign_int(ctx->K, &pp->value.num, 1);
      }
    }
    break;
  case COEFFICIENT_POLYNOMIAL:
  {
    int i;
    coefficient_t gcd;
    // Compute the gcd of coefficients starting with LC
    coefficient_construct_copy(ctx, &gcd, coefficient_lc(C));
    // Make if positive in case it's the only one
    if (coefficient_lc_sgn(ctx, &gcd) < 0) {
      coefficient_neg(ctx, &gcd, &gcd);
    }
    // Compute the rest of the gcd
    for (i = SIZE(C)-2; i >= 0 ; -- i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_gcd(ctx, &gcd, &gcd, COEFF(C, i));
        if (coefficient_is_one(ctx, &gcd)) {
          break;
        }
      }
    }
    // GCD is positive, so if the leading coefficient of C is negative, flip it
    if (coefficient_lc_sgn(ctx, C) < 0) {
      coefficient_neg(ctx, &gcd, &gcd);
    }

    if (pp) {
      // Now compute the pp
      coefficient_div(ctx, pp, C, &gcd);
      assert(coefficient_is_normalized(ctx, pp));
    }
    if (cont) {
      coefficient_swap(&gcd, cont);
      assert(coefficient_is_normalized(ctx, cont));
    }
    coefficient_destruct(&gcd);
    break;
  }
  default:
    assert(0);
    break;
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_pp_cont() => ");
    if (pp) { tracef("pp = "); coefficient_ops.print(ctx, pp, trace_out); tracef("\n"); }
    if (cont) { tracef("cont = "); coefficient_ops.print(ctx, cont, trace_out); tracef("\n"); }
  }
}

static
void coefficient_cont(const polynomial_context_t* ctx, coefficient_t* cont, const coefficient_t* C) {
  coefficient_pp_cont(ctx, 0, cont, C);
}

static
void coefficient_pp(const polynomial_context_t* ctx, coefficient_t* pp, const coefficient_t* C) {
  coefficient_pp_cont(ctx, pp, 0, C);
}


void monomial_gcd_visit(const polynomial_context_t* ctx, monomial_t* m, void* data) {
  monomial_t* gcd = (monomial_t*) data;
  if (integer_is_zero(ctx->K, &gcd->a)) {
    monomial_ops.assign(ctx, gcd, m, 0);
  } else {
    monomial_ops.gcd(ctx, gcd, gcd, m);
  }
}

static
int coefficient_is_univariate(const polynomial_context_t* ctx, const coefficient_t* C) {
  int i;
  if (C->type == COEFFICIENT_NUMERIC) {
    return 1;
  } else {
    for (i = 0; i < SIZE(C); ++i) {
      if (COEFF(C, i)->type != COEFFICIENT_NUMERIC) {
        return 0;
      }
    }
    return 1;
  }
}

static
const integer_t* coefficient_get_constant(const coefficient_t* C) {
  while (C->type == COEFFICIENT_POLYNOMIAL) {
    C = COEFF(C, 0);
  }
  return &C->value.num;
}

static
upolynomial_t* coefficient_to_simple_univariate(const polynomial_context_t* ctx, const coefficient_t* C) {
  assert(C->type == COEFFICIENT_POLYNOMIAL);

  integer_t* coeff = malloc(sizeof(integer_t)*SIZE(C));

  int i;
  for (i = 0; i < SIZE(C); ++ i) {
    integer_construct_copy(ctx->K, coeff + i, coefficient_get_constant(COEFF(C, i)));
  }

  upolynomial_t* C_u = upolynomial_ops.construct(ctx->K, SIZE(C) - 1, coeff);

  for (i = 0; i < SIZE(C); ++ i) {
    integer_destruct(coeff + i);
  }

  return C_u;
}

/**
 * Takes two (primitive) coefficients over the same variable, makes them univariate by
 * substituting 0 for other variables (if any). Then it computes the
 * univariate gcd of these. If the coefficients were univariate already, or
 * the result is a constant (i.e. gcd = 1), the result is precise.
 */
static
int coefficient_gcd_pp_univariate(const polynomial_context_t* ctx,
    coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2) {

  assert(C1->type == COEFFICIENT_POLYNOMIAL);
  assert(C2->type == COEFFICIENT_POLYNOMIAL);

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp_univariate()\n");
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  int C1_vanishes = integer_is_zero(ctx->K, coefficient_get_constant(coefficient_lc(C1)));
  int C2_vanishes = integer_is_zero(ctx->K, coefficient_get_constant(coefficient_lc(C2)));

  if (C1_vanishes || C2_vanishes) {
    // One of C1 or C2 vanishes in the univariate conversion, we're not precise enough
    return 0;
  }

  variable_t x = VAR(C1);
  assert(x == VAR(C2));

  upolynomial_t* C1_u = coefficient_to_simple_univariate(ctx, C1);
  upolynomial_t* C2_u = coefficient_to_simple_univariate(ctx, C2);
  upolynomial_t* gcd_u = upolynomial_ops.gcd(C1_u, C2_u);
  coefficient_construct_from_univariate(ctx, gcd, gcd_u, x);

  upolynomial_ops.destruct(C1_u);
  upolynomial_ops.destruct(C2_u);
  upolynomial_ops.destruct(gcd_u);

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp_univariate() => ");
    tracef("gcd = "); coefficient_ops.print(ctx, gcd, trace_out); tracef("\n");
  }

  if (gcd->type == COEFFICIENT_NUMERIC) {
    integer_assign_int(ctx->K, &gcd->value.num, 1);
    return 1;
  } else if (coefficient_is_univariate(ctx, C1) && coefficient_is_univariate(ctx, C2)) {
    return 1;
  } else {
    return 0;
  }
}

/**
 * Extracts the largest monomial power out of P and Q and into gcd, also divide.
 * For example, P and Q in Z[y, x]
 *
 *  P = 4*y*x^2 + 2*y^2 = 2*y^2*(2*x^2 + 1)
 *  Q = 2*y^3*x^3
 *
 * gives
 *
 *  gcd = 2*y^2
 *  P = 2*x^2 + 1
 *  Q = 2*y*x^3
 */
static
void coefficient_gcd_monomial_extract(const polynomial_context_t* ctx, coefficient_t* gcd, coefficient_t* P, coefficient_t* Q) {

  TRACE("coefficient", "coefficient_gcd_monomial_extract()\n");

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_ops.print(ctx, Q, trace_out); tracef("\n");
  }

  assert(P != Q);

  monomial_t m_P_gcd, m_Q_gcd, m_tmp;
  monomial_ops.construct(ctx, &m_P_gcd);
  monomial_ops.construct(ctx, &m_Q_gcd);
  monomial_ops.construct(ctx, &m_tmp);

  // Compute the gcd
  coefficient_traverse(ctx, P, monomial_gcd_visit, &m_tmp, &m_P_gcd);
  monomial_ops.clear(ctx, &m_tmp);
  coefficient_traverse(ctx, Q, monomial_gcd_visit, &m_tmp, &m_Q_gcd);

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("P_gcd = "); monomial_ops.print(ctx, &m_P_gcd, trace_out); tracef("\n");
    tracef("Q_gcd = "); monomial_ops.print(ctx, &m_Q_gcd, trace_out); tracef("\n");
  }

  // Final gcd
  monomial_t m_gcd;
  monomial_ops.construct(ctx, &m_gcd);
  monomial_ops.gcd(ctx, &m_gcd, &m_P_gcd, &m_Q_gcd);

  // Construct the result
  coefficient_t result;
  coefficient_construct(ctx, &result);
  coefficient_add_monomial(ctx, &m_gcd, &result);

  // Divide P and Q with their gcds
  coefficient_t P_gcd, Q_gcd;
  coefficient_construct(ctx, &P_gcd);
  coefficient_construct(ctx, &Q_gcd);
  coefficient_add_monomial(ctx, &m_P_gcd, &P_gcd);
  coefficient_add_monomial(ctx, &m_Q_gcd, &Q_gcd);
  coefficient_div(ctx, P, P, &P_gcd);
  coefficient_div(ctx, Q, Q, &Q_gcd);
  coefficient_destruct(&P_gcd);
  coefficient_destruct(&Q_gcd);

  // Output the result
  coefficient_swap(&result, gcd);
  coefficient_destruct(&result);

  monomial_ops.destruct(&m_gcd);
  monomial_ops.destruct(&m_tmp);
  monomial_ops.destruct(&m_Q_gcd);
  monomial_ops.destruct(&m_P_gcd);

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_gcd_monomial_extract() =>"); coefficient_ops.print(ctx, gcd, trace_out); tracef("\n");
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_ops.print(ctx, Q, trace_out); tracef("\n");
  }
}

STAT_DECLARE(int, coefficient, gcd_pp);

/**
 * Compute the gcd of two primitive polynomials P and Q. The polynomials P and
 * Q will be used and changed in the computation.
 */
static
void coefficient_gcd_pp(const polynomial_context_t* ctx, coefficient_t* gcd, coefficient_t* P, coefficient_t* Q) {

  TRACE("coefficient", "coefficient_gcd_pp()\n");
  STAT(coefficient, gcd_pp) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_ops.print(ctx, Q, trace_out); tracef("\n");
  }

  // Try to comute the univariate GCD first
  coefficient_t gcd_u;
  coefficient_construct(ctx, &gcd_u);

  int precise = coefficient_gcd_pp_univariate(ctx, &gcd_u, P, Q);
  if (precise) {
    // GCD = 1, just copy the univariate gcd
    coefficient_swap(gcd, &gcd_u);
  } else {

    coefficient_t R;
    coefficient_construct(ctx, &R);

    //
    // We compute the reduction of P and Q in Z[y, x], i.e.
    //
    //   a*P = b*Q + R
    //
    // with a in Z[y], b in Z[y, x], and deg(R) < deg(Q) or deg(R) == 0.
    //
    // P and Q are primitive so GCD(P, Q) should be primitive, i.e. in
    // Z[y, x] or 1. therefore GCD(P, Q) = 1 if R != 0, or
    // GCD(P, Q) = po(Q) if R = 0
    //
    do {

      // One step reduction
      coefficient_reduce(ctx, P, Q, 0, 0, &R, REMAINDERING_LCM_SPARSE);

      int cmp_type = coefficient_cmp_type(ctx, Q, &R);
      if (cmp_type == 0) {
        // P = Q
        // Q = pp(R)
        coefficient_swap(P, Q);
        coefficient_swap(Q, &R);
        coefficient_pp(ctx, Q, Q);
      } else {
        assert(cmp_type > 0);
        if (!coefficient_is_zero(ctx, &R)) {
          coefficient_destruct(Q);
          coefficient_construct_from_int(ctx, Q, 1);
        }
        break;
      }
    } while (1);

    coefficient_swap(Q, gcd);
    coefficient_destruct(&R);
  }

  coefficient_destruct(&gcd_u);

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_gcd_pp() => "); coefficient_ops.print(ctx, gcd, trace_out); tracef("\n");
  }
}

STAT_DECLARE(int, coefficient, gcd);

static
void coefficient_gcd(const polynomial_context_t* ctx, coefficient_t* gcd, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_gcd()\n");
  STAT(coefficient, gcd) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(ctx->K == Z);

  int cmp_type = coefficient_cmp_type(ctx, C1, C2);

  if (cmp_type < 0) {
    const coefficient_t* tmp = C1;
    C1 = C2;
    C2 = tmp;
    cmp_type = -cmp_type;
  }

  if (cmp_type == 0) {
    switch (C1->type) {
    case COEFFICIENT_NUMERIC:
      if (gcd->type == COEFFICIENT_POLYNOMIAL) {
        coefficient_destruct(gcd);
        coefficient_construct(ctx, gcd);
      }
      integer_gcd_Z(&gcd->value.num, &C1->value.num, &C2->value.num);
      break;
    case COEFFICIENT_POLYNOMIAL:
    {
      coefficient_t P, Q;
      if (SIZE(C1) > SIZE(C2)) {
        coefficient_construct_copy(ctx, &P, C1);
        coefficient_construct_copy(ctx, &Q, C2);
      } else {
        coefficient_construct_copy(ctx, &P, C2);
        coefficient_construct_copy(ctx, &Q, C1);
      }

      // Get the common power variables out
      coefficient_t gcd_mon;
      coefficient_construct(ctx, &gcd_mon);
      coefficient_gcd_monomial_extract(ctx, &gcd_mon, &P, &Q);

      // If monomial extraction changed the type, we need to go again
      if (coefficient_cmp_type(ctx, C1, &P) != 0 || coefficient_cmp_type(ctx, C2, &Q) != 0) {
        coefficient_gcd(ctx, gcd, &P, &Q);
      } else {
        // Normalize the P and Q to be primitive (and keep the content)
        coefficient_t P_cont, Q_cont;
        coefficient_construct(ctx, &P_cont);
        coefficient_construct(ctx, &Q_cont);
        coefficient_pp_cont(ctx, &P, &P_cont, &P);
        coefficient_pp_cont(ctx, &Q, &Q_cont, &Q);

        // Get the gcd of the content
        coefficient_t gcd_cont;
        coefficient_construct(ctx, &gcd_cont);
        coefficient_gcd(ctx, &gcd_cont, &P_cont, &Q_cont);

        // Get the gcd of the primitive parts
        coefficient_gcd_pp(ctx, gcd, &P, &Q);

        // Multiply in the content gcd
        coefficient_mul(ctx, gcd, gcd, &gcd_cont);

        coefficient_destruct(&P_cont);
        coefficient_destruct(&Q_cont);
        coefficient_destruct(&gcd_cont);
      }

      // Multiply in the monomial gcd
      coefficient_mul(ctx, gcd, gcd, &gcd_mon);

      // Remove temps
      coefficient_destruct(&P);
      coefficient_destruct(&Q);
      coefficient_destruct(&gcd_mon);
      break;
    }
    default:
      assert(0);
      break;
    }
  } else {
    // C1 in Z[y, x]
    // C2 in Z[y]
    // so GCD(C1, C2) = GCD(cont(C1), C2)
    coefficient_t cont;
    coefficient_construct(ctx, &cont);
    coefficient_cont(ctx, &cont, C1);
    coefficient_gcd(ctx, gcd, &cont, C2);
    coefficient_destruct(&cont);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_gcd() => "); coefficient_ops.print(ctx, gcd, trace_out); tracef("\n");
  }

  if (debug_trace_ops.is_enabled("coefficient::gcd::sage")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
    tracef("gcd = "); coefficient_ops.print(ctx, gcd, trace_out); tracef("\n");
    tracef("gcd_sage = C1.gcd(C2)\n");
    tracef("if (gcd != gcd_sage):\n");
    tracef("\tprint 'C1 =', C1\n");
    tracef("\tprint 'C2 =', C2\n");
  }

  assert(coefficient_is_normalized(ctx, gcd));
}

STAT_DECLARE(int, coefficient, lcm);

static
void coefficient_lcm(const polynomial_context_t* ctx, coefficient_t* lcm, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_lcm()\n");
  STAT(coefficient, lcm) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_ops.print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_ops.print(ctx, C2, trace_out); tracef("\n");
  }

  assert(ctx->K == Z);

  if (C1->type == COEFFICIENT_NUMERIC && C2->type == COEFFICIENT_NUMERIC) {
    // Integer LCM
    if (lcm->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_destruct(lcm);
      coefficient_construct(ctx, lcm);
    }
    integer_lcm_Z(&lcm->value.num, &C1->value.num, &C2->value.num);
  } else {
    // LCM(C1, C2) = C1*C2/GCD(C1, C2)
    coefficient_t gcd;
    coefficient_construct(ctx, &gcd);
    coefficient_gcd(ctx, &gcd, C1, C2);
    if (coefficient_is_one(ctx, &gcd)) {
      coefficient_mul(ctx, lcm, C1, C2);
    } else {
      if (coefficient_cmp_type(ctx, C1, C2) <= 0) {
        coefficient_div(ctx, lcm, C1, &gcd);
        coefficient_mul(ctx, lcm, lcm, C2);
      } else {
        coefficient_div(ctx, lcm, C2, &gcd);
        coefficient_mul(ctx, lcm, lcm, C1);
      }
    }
    if (coefficient_lc_sgn(ctx, lcm) < 0) {
      coefficient_neg(ctx, lcm, lcm);
    }
    coefficient_destruct(&gcd);
  }

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("coefficient_lcm() => "); coefficient_ops.print(ctx, lcm, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, lcm));
}

STAT_DECLARE(int, coefficient, psc);

/**
 * (non-optimized) Subresultant algorithm, as described in
 *
 * [2000] Ducos - Optimizations of the subresltant algorithm.
 */
void coefficient_psc(const polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* P, const coefficient_t* Q) {

  TRACE("coefficient", "coefficient_psc()\n");
  STAT(coefficient, psc) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("P = "); coefficient_ops.print(ctx, P, trace_out); tracef("\n");
    tracef("Q = "); coefficient_ops.print(ctx, Q, trace_out); tracef("\n");
  }

  assert(P->type == COEFFICIENT_POLYNOMIAL);
  assert(Q->type == COEFFICIENT_POLYNOMIAL);

  variable_t x = VAR(P);
  assert(VAR(Q) == x);

  size_t P_deg = coefficient_degree(P);
  size_t Q_deg = coefficient_degree(Q);
  assert(P_deg >= Q_deg);

  // S = []
  int S_size = 0;

  // s = lc(Q)^(deg(P) - deg(Q)
  coefficient_t s;
  coefficient_construct(ctx, &s);
  coefficient_pow(ctx, &s, coefficient_lc(Q), P_deg - Q_deg);

  // Set the final position
  coefficient_assign(ctx, S + Q_deg, &s);

  // A = Q, B = prem(P, -Q)
  coefficient_t A, B;
  coefficient_construct_copy(ctx, &A, Q);
  coefficient_construct_copy(ctx, &B, Q);
  coefficient_neg(ctx, &B, &B);
  coefficient_prem(ctx, &B, P, &B);

  // Some temporaries
  coefficient_t C, pow;
  coefficient_construct(ctx, &C);
  coefficient_construct(ctx, &pow);

  for (;;) {

    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
    }

    // d = deg(A); e = deg(B)
    size_t d = coefficient_degree_safe(ctx, &A, x);
    size_t e = coefficient_degree_safe(ctx, &B, x);
    assert(d > e);

    // Holds:
    //   A ~ S_d   if d = deg(Q)
    //   A = S_d   if d < deg(Q)
    //   B = S_d-1, s = lc(S_d) for d <= deg(Q)

    if (coefficient_is_zero(ctx, &B)) {
      break;
    }

    // S = [B; S]
    if (coefficient_degree_safe(ctx, &B, x) == d - 1) {
      coefficient_assign(ctx, S + S_size, coefficient_lc_safe(ctx, &B, x));
    }
    S_size ++;
    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("S[%d] = ", S_size - 1); coefficient_print(ctx, S + S_size - 1, trace_out); tracef("\n");
    }

    // Holds:
    //   S = [S_d-1, S_d, ...]

    int delta = d - e;
    if (delta > 1) {
      // C = (lc(B)^(delta-1)*B)/(s^(delta-1))
      coefficient_pow(ctx, &pow, coefficient_lc_safe(ctx, &B, x), delta-1);
      coefficient_mul(ctx, &C, &pow, &B);
      coefficient_pow(ctx, &pow, &s, delta-1);
      coefficient_div(ctx, &C, &C, &pow);
      // S = [C; S]
      if (coefficient_degree_safe(ctx, &C, x) == e) {
        coefficient_assign(ctx, S + S_size, coefficient_lc_safe(ctx, &C, x));
      }
      S_size ++;
      if (debug_trace_ops.is_enabled("coefficient::resultant")) {
        tracef("S[%d] = ", S_size - 1); coefficient_print(ctx, S + S_size - 1, trace_out); tracef("\n");
      }
    } else {
      // C = B
      coefficient_assign(ctx, &C, &B);
    }

    // Holds:
    //   C = S_e, S = [S_e, ...]

    if (e == 0) {
      break;
    }

    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("lc(A) = "); coefficient_print(ctx, coefficient_lc_safe(ctx, &A, x), trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
      tracef("s = "); coefficient_print(ctx, &s, trace_out); tracef("\n");
      tracef("delta = %d\n", delta);
    }

    // B = prem(A, -B)/(s^delta*lc(A))
    coefficient_neg(ctx, &B, &B);
    coefficient_prem(ctx, &B, &A, &B);
    coefficient_pow(ctx, &pow, &s, delta);
    coefficient_mul(ctx, &pow, &pow, coefficient_lc_safe(ctx, &A, x));
    coefficient_div(ctx, &B, &B, &pow);

    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("B = "); coefficient_print(ctx, &B, trace_out); tracef("\n");
    }

    // Holds:
    //   B = S_e-1

    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("C = "); coefficient_print(ctx, &C, trace_out); tracef("\n");
      tracef("lc(A) = "); coefficient_print(ctx, coefficient_lc_safe(ctx, &A, x), trace_out); tracef("\n");
    }

    // A = C, s = lc(A)
    coefficient_swap(&A, &C);
    coefficient_assign(ctx, &s, coefficient_lc_safe(ctx, &A, x));

    if (debug_trace_ops.is_enabled("coefficient::resultant")) {
      tracef("A = "); coefficient_print(ctx, &A, trace_out); tracef("\n");
      tracef("s = "); coefficient_print(ctx, &s, trace_out); tracef("\n");
    }
  }

  // Reverse S
  int i = 0, j = S_size - 1;
  while (i < j) {
    coefficient_swap(S + i, S + j);
    i ++;
    j --;
  }

  // Remove temps
  coefficient_destruct(&A);
  coefficient_destruct(&B);
  coefficient_destruct(&C);
  coefficient_destruct(&pow);
  coefficient_destruct(&s);
}

STAT_DECLARE(int, coefficient, resultant);

void coefficient_resultant(const polynomial_context_t* ctx, coefficient_t* res, const coefficient_t* A, const coefficient_t* B) {
  TRACE("coefficient", "coefficient_resultant()\n");
  STAT(coefficient, resultant) ++;

  if (debug_trace_ops.is_enabled("coefficient")) {
    tracef("A = "); coefficient_ops.print(ctx, A, trace_out); tracef("\n");
    tracef("B = "); coefficient_ops.print(ctx, B, trace_out); tracef("\n");
  }

  assert(A->type == COEFFICIENT_POLYNOMIAL);
  assert(B->type == COEFFICIENT_POLYNOMIAL);

  assert(VAR(B) == VAR(A));

  size_t A_deg = coefficient_degree(A);
  size_t B_deg = coefficient_degree(B);

  if (A_deg < B_deg) {
    coefficient_resultant(ctx, res, B, A);
    if ((A_deg % 2) && (B_deg % 2)) {
      coefficient_neg(ctx, res, res);
    }
    return;
  }

  // Compute the PSC
  size_t psc_size = B_deg + 1;
  coefficient_t* psc = malloc(sizeof(coefficient_t)*psc_size);
  int i;
  for (i = 0; i < psc_size; ++ i) {
    coefficient_construct(ctx, psc + i);
  }
  coefficient_psc(ctx, psc, A, B);

  // Resultant is PSC[0]
  coefficient_swap(res, psc);

  // Remove temps
  for (i = 0; i < psc_size; ++ i) {
    coefficient_destruct(psc + i);
  }
  free(psc);
}

///
/// Normalization that everyone is using
///

static void
coefficient_normalize(const polynomial_context_t* ctx, coefficient_t* C) {
  if (C->type == COEFFICIENT_POLYNOMIAL) {
    assert(C->value.rec.size >= 1);
    size_t i = C->value.rec.size - 1;
    // Find the first non-zero coefficient
    while (i > 0 && coefficient_is_zero(ctx, COEFF(C, i))) {
      i --;
    }
    // If a constant, just upgrade it (it could be an actual constant or just
    // some other polynomial)
    if (i == 0) {
      coefficient_t result;
      coefficient_construct(ctx, &result);
      coefficient_swap(&result, COEFF(C, 0));
      coefficient_swap(&result, C);
      coefficient_destruct(&result);
    } else {
      // We are a proper polynomial, set the size
      C->value.rec.size = i + 1;
    }
  }
}

int
coefficient_is_normalized(const polynomial_context_t* ctx, coefficient_t* C) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 1;
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (SIZE(C) <= 1) {
      return 0;
    }
    if (coefficient_is_zero(ctx, COEFF(C, SIZE(C) - 1))) {
      return 0;
    }
    return 1;
    break;
  default:
    assert(0);
    break;
  }
  return 0;
}

static void
coefficient_ensure_capacity(const polynomial_context_t* ctx, coefficient_t* C, variable_t x, size_t capacity) {
  assert(capacity >= 1);
  coefficient_t tmp;
  switch(C->type) {
  case COEFFICIENT_NUMERIC:
    // Create a new recursive and put the current constant into constants place
    coefficient_construct_rec(ctx, &tmp, x, capacity);
    coefficient_swap(COEFF(&tmp, 0), C);
    coefficient_swap(C, &tmp);
    coefficient_destruct(&tmp);
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (x != VAR(C)) {
      assert(ctx->var_order->ops->cmp(ctx->var_order, x, VAR(C)) > 0);
      // Same as for constants above
      coefficient_construct_rec(ctx, &tmp, x, capacity);
      coefficient_swap(COEFF(&tmp, 0), C);
      coefficient_swap(C, &tmp);
      coefficient_destruct(&tmp);
    } else if (capacity > C->value.rec.capacity) {
      // Already recursive polynomial, resize up
      C->value.rec.coefficients = realloc(C->value.rec.coefficients, capacity * sizeof(coefficient_t));
      int i;
      for (i = C->value.rec.capacity; i < capacity; ++ i) {
        coefficient_construct(ctx, C->value.rec.coefficients + i);
      }
      C->value.rec.capacity = capacity;
      C->value.rec.size = capacity;
    }
    break;
  }
}



const coefficient_ops_t coefficient_ops = {
    coefficient_construct,
    coefficient_construct_from_int,
    coefficient_construct_from_integer,
    coefficient_construct_from_univariate,
    coefficient_construct_simple,
    coefficient_construct_copy,
    coefficient_destruct,
    coefficient_swap,
    coefficient_assign,
    coefficient_assign_int,
    coefficient_to_simple_univariate,
    coefficient_is_constant,
    coefficient_degree,
    coefficient_top_variable,
    coefficient_get_coefficient,
    coefficient_lc,
    coefficient_is_zero,
    coefficient_is_one,
    coefficient_sgn,
    coefficient_value_approx,
    coefficient_lc_sgn,
    coefficient_in_order,
    coefficient_cmp,
    coefficient_cmp_type,
    coefficient_divides,
    coefficient_print,
    coefficient_to_string,
    coefficient_order,
    coefficient_add,
    coefficient_sub,
    coefficient_neg,
    coefficient_mul,
    coefficient_mul_int,
    coefficient_shl,
    coefficient_shr,
    coefficient_pow,
    coefficient_add_mul,
    coefficient_sub_mul,
    coefficient_reduce,
    coefficient_div,
    coefficient_rem,
    coefficient_prem,
    coefficient_sprem,
    coefficient_divrem,
    coefficient_derivative,
    coefficient_gcd,
    coefficient_pp,
    coefficient_cont,
    coefficient_pp_cont,
    coefficient_lcm,
    coefficient_resultant,
    coefficient_psc,
    coefficient_set_power_symbol
};
