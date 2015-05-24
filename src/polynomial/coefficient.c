/*
 * polynomial_coefficient.c
 *
 *  Created on: Feb 4, 2014
 *      Author: dejan
 */

#include <interval.h>
#include <upolynomial.h>
#include <variable_db.h>
#include <monomial.h>
#include <variable_list.h>

#include "polynomial/polynomial.h"
#include "polynomial/coefficient.h"
#include "polynomial/output.h"
#include "polynomial/gcd.h"

#include "upolynomial/upolynomial.h"

#include "number/integer.h"
#include "number/rational.h"
#include "interval/arithmetic.h"

#include "variable/variable_order.h"
#include "polynomial/polynomial_context.h"

#include <assignment.h>

#include "utils/debug_trace.h"
#include "utils/statistics.h"

#include <assert.h>

static
void coefficient_resolve_algebraic(const lp_polynomial_context_t* ctx, const coefficient_t* A, const lp_assignment_t* m, coefficient_t* A_alg);

static
lp_dyadic_interval_t* algebraic_interval_remember(const lp_variable_list_t* var_list, const lp_assignment_t* m) {
  // Remember the existing intervals so that we can recover later
  lp_dyadic_interval_t* x_i_intervals = malloc(sizeof(lp_dyadic_interval_t)*var_list->list_size);
  size_t i;
  for (i = 0; i < var_list->list_size; ++ i) {
    lp_variable_t x_i = var_list->list[i];
    const lp_value_t* x_i_value = lp_assignment_get_value(m, x_i);
    if (!lp_value_is_rational(x_i_value)) {
      lp_dyadic_interval_construct_copy(x_i_intervals + i, &x_i_value->value.a.I);
    } else {
      lp_dyadic_interval_construct_zero(x_i_intervals + i);
    }
  }
  return x_i_intervals;
}

static
void algebraic_interval_restore(const lp_variable_list_t* var_list, lp_dyadic_interval_t* cache, const lp_assignment_t* m) {
  // Restore the intervals, and destroy the temps
  size_t i;
  for (i = 0; i < var_list->list_size; ++ i) {
    lp_variable_t x_i = var_list->list[i];
    const lp_value_t* x_i_value = lp_assignment_get_value(m, x_i);
    if (x_i_value->type == LP_VALUE_ALGEBRAIC && !x_i_value->value.a.I.is_point) {
      lp_algebraic_number_restore_interval_const(&x_i_value->value.a, cache + i);
    }
    lp_dyadic_interval_destruct(cache + i);
  }
  free(cache);
}


/**
 * Make sure that the coefficient has the given capacity for the given variable.
 * If the polynomial a constant or a polynomial in a smaller variable, we
 * transform it into a polynomial and adjust the coefficent.
 */
static void
coefficient_ensure_capacity(const lp_polynomial_context_t* ctx, coefficient_t* C, lp_variable_t x, size_t capacity);

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
coefficient_normalize(const lp_polynomial_context_t* ctx, coefficient_t* C);

int
coefficient_is_normalized(const lp_polynomial_context_t* ctx, coefficient_t* C);

STAT_DECLARE(int, coefficient, construct)

void coefficient_construct(const lp_polynomial_context_t* ctx, coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_construct()\n");
  STAT(coefficient, construct) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_from_int(ctx->K, &C->value.num, 0);
}

STAT_DECLARE(int, coefficient, construct_from_int)

void coefficient_construct_from_int(const lp_polynomial_context_t* ctx, coefficient_t* C, long C_int) {
  TRACE("coefficient::internal", "coefficient_construct_from_int()\n");
  STAT(coefficient, construct_from_int) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_from_int(ctx->K, &C->value.num, C_int);
}

STAT_DECLARE(int, coefficient, construct_from_integer)

void coefficient_construct_from_integer(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* C_integer) {
  TRACE("coefficient::internal", "coefficient_construct_from_integer()\n");
  STAT(coefficient, construct_from_integer) ++;

  C->type = COEFFICIENT_NUMERIC;
  integer_construct_copy(ctx->K, &C->value.num, C_integer);
}

STAT_DECLARE(int, coefficient, construct_rec)

void coefficient_construct_rec(const lp_polynomial_context_t* ctx, coefficient_t* C, lp_variable_t x, size_t capacity) {
  TRACE("coefficient::internal", "coefficient_construct_rec()\n");
  STAT(coefficient, construct_rec) ++;

  C->type = COEFFICIENT_POLYNOMIAL;
  C->value.rec.x = x;
  C->value.rec.size = 0;
  C->value.rec.capacity = 0;
  C->value.rec.coefficients = 0;
  coefficient_ensure_capacity(ctx, C, x, capacity);
}

STAT_DECLARE(int, coefficient, construct_simple)

void coefficient_construct_simple(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_integer_t* a, lp_variable_t x, unsigned n) {
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

STAT_DECLARE(int, coefficient, construct_copy)

void coefficient_construct_copy(const lp_polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from) {
  TRACE("coefficient::internal", "coefficient_construct_copy()\n");
  STAT(coefficient, construct_copy) ++;

  size_t i;
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

STAT_DECLARE(int, coefficient, construct_from_univariate)

void coefficient_construct_from_univariate(const lp_polynomial_context_t* ctx,
    coefficient_t* C, const lp_upolynomial_t* C_u, lp_variable_t x) {

  TRACE("coefficient::internal", "coefficient_construct_from_univariate()\n");
  STAT(coefficient, construct_from_univariate) ++;

  // Get the coefficients
  size_t C_u_deg = lp_upolynomial_degree(C_u);
  lp_integer_t* coeff = malloc(sizeof(lp_integer_t)*(C_u_deg + 1));

  size_t i;
  for (i = 0; i <= C_u_deg; ++ i) {
    integer_construct_from_int(ctx->K, coeff + i, 0);
  }

  lp_upolynomial_unpack(C_u, coeff);

  // Construct the polynomial
  coefficient_construct_rec(ctx, C, x, C_u_deg + 1);

  // Move over the coefficients
  for (i = 0; i <= C_u_deg; ++ i) {
    integer_swap(&COEFF(C, i)->value.num, coeff + i);
    integer_destruct(coeff + i);
  }
  free(coeff);

  // Normalize (it might be constant)
  coefficient_normalize(ctx, C);

  assert(coefficient_is_normalized(ctx, C));
}

void coefficient_destruct(coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_destruct()\n");

  size_t i;
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    integer_destruct(&C->value.num);
    break;
  case COEFFICIENT_POLYNOMIAL:
    for (i = 0; i < CAPACITY(C); ++ i) {
      coefficient_destruct(COEFF(C, i));
    }
    free(C->value.rec.coefficients);
    break;
  default:
    assert(0);
  }
}

STAT_DECLARE(int, coefficient, swap)

void coefficient_swap(coefficient_t* C1, coefficient_t* C2) {
  TRACE("coefficient::internal", "coefficient_swap()\n");
  STAT(coefficient, swap) ++;
  coefficient_t tmp = *C1;
  *C1 = *C2;
  *C2 = tmp;
}

STAT_DECLARE(int, coefficient, assign)

void coefficient_assign(const lp_polynomial_context_t* ctx, coefficient_t* C, const coefficient_t* from) {
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

STAT_DECLARE(int, coefficient, assign_int)

void coefficient_assign_int(const lp_polynomial_context_t* ctx, coefficient_t* C, long x) {
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

const coefficient_t* coefficient_lc_safe(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_variable_t x) {
  __unused(ctx);
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return C;
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (VAR(C) == x) {
      return COEFF(C, SIZE(C) - 1);
    } else {
      assert(lp_variable_order_cmp(ctx->var_order, x, VAR(C)) > 0);
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

const coefficient_t* coefficient_lc_m(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return C;
    break;
  case COEFFICIENT_POLYNOMIAL: {
    // Locate the first non-zero coefficient past the top one
    int i = SIZE(C) - 1;
    while (i > 0 && coefficient_sgn(ctx, COEFF(C, i), M) == 0) {
      -- i;
    }
    return COEFF(C, i);
    break;
  }
  }
  assert(0);
  return 0;
}


void coefficient_reductum(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C) {

  assert(C->type == COEFFICIENT_POLYNOMIAL);

  // Locate the first non-zero ceofficient past the top one
  int i = SIZE(C) - 2;
  while (coefficient_is_zero(ctx, COEFF(C, i))) {
    -- i;
  }

  coefficient_t result;
  coefficient_construct_rec(ctx, &result, VAR(C), i + 1);

  // Copy the other coefficients
  while (i >= 0) {
    if (!coefficient_is_zero(ctx, COEFF(C, i))) {
      coefficient_assign(ctx, COEFF(&result, i), COEFF(C, i));
    }
    -- i;
  }

  coefficient_normalize(ctx, &result);
  coefficient_swap(R, &result);
  coefficient_destruct(&result);
}

void coefficient_reductum_m(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C, const lp_assignment_t* m) {

  assert(C->type == COEFFICIENT_POLYNOMIAL);

  // Locate the first non-zero ceofficient past the top one
  int i = SIZE(C) - 2;
  while (i > 0 && coefficient_sgn(ctx, COEFF(C, i), m) == 0) {
    -- i;
  }

  coefficient_t result;
  coefficient_construct_rec(ctx, &result, VAR(C), i + 1);

  // Copy the other coefficients
  while (i >= 0) {
    if (!coefficient_is_zero(ctx, COEFF(C, i))) {
      coefficient_assign(ctx, COEFF(&result, i), COEFF(C, i));
    }
    -- i;
  }

  coefficient_normalize(ctx, &result);
  coefficient_swap(R, &result);
  coefficient_destruct(&result);
}

int coefficient_is_constant(const coefficient_t* C) {
  return C->type == COEFFICIENT_NUMERIC;
}

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

size_t coefficient_degree_m(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M) {
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 0;
    break;
  case COEFFICIENT_POLYNOMIAL: {
    // Locate the first non-zero coefficient past the top one
    size_t i = SIZE(C) - 1;
    while (i > 0 && coefficient_sgn(ctx, COEFF(C, i), M) == 0) {
      -- i;
    }
    return i;
    break;
  }
  }
  assert(0);
  return 0;

}

size_t coefficient_degree_safe(const lp_polynomial_context_t* ctx, const coefficient_t* C, lp_variable_t x) {
  __unused(ctx);
  switch (C->type) {
  case COEFFICIENT_NUMERIC:
    return 0;
    break;
  case COEFFICIENT_POLYNOMIAL:
    if (VAR(C) == x) {
      return SIZE(C) - 1;
    } else {
      assert(lp_variable_order_cmp(ctx->var_order, x, VAR(C)) > 0);
      return 0;
    }
    break;
  }
  assert(0);
  return 0;
}

lp_variable_t coefficient_top_variable(const coefficient_t* C) {
  if (C->type == COEFFICIENT_POLYNOMIAL) {
    return VAR(C);
  } else {
    return lp_variable_null;
  }
}

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

STAT_DECLARE(int, coefficient, is_zero)

int coefficient_is_zero(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, is_zero) ++;
  return C->type == COEFFICIENT_NUMERIC && integer_is_zero(ctx->K, &C->value.num);
}

STAT_DECLARE(int, coefficient, is_one)

int coefficient_is_one(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, is_one) ++;
  return C->type == COEFFICIENT_NUMERIC && integer_cmp_int(ctx->K, &C->value.num, 1) == 0;
}

void coefficient_value_approx(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m, lp_rational_interval_t* value) {

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_value_approx("); coefficient_print(ctx, C, trace_out); tracef(")\n");
  }

  if (C->type == COEFFICIENT_NUMERIC) {
    lp_rational_interval_t result;
    lp_rational_interval_construct_from_integer(&result, &C->value.num, 0, &C->value.num, 0);
    lp_rational_interval_swap(value, &result);
    lp_rational_interval_destruct(&result);
  } else {

    lp_rational_interval_t result, tmp1, tmp2, x_value;

    lp_rational_interval_construct_zero(&result);
    lp_rational_interval_construct_zero(&tmp1);
    lp_rational_interval_construct_zero(&tmp2);
    lp_rational_interval_construct_zero(&x_value);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_value_approx(): x = %s\n", lp_variable_db_get_name(ctx->var_db, VAR(C)));
    }

    lp_assignment_get_value_approx(m, VAR(C), &x_value);

    // Get the value of x
    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_value_approx(): x_value = ");
      lp_rational_interval_print(&x_value, trace_out);
      tracef("\n");
    }

    // We compute using powers, just an attempt to compute better. For example
    // if p = x^2 + x and x = [-1, 1] then
    //  a) x(x + 1) = [-1, 1]*[0, 2] = [-2, 2]
    //  b) x^2 + x = [0, 1] + [-1, 1] = [-1, 2]
    // we choose to do it b) way

    // Compute
    size_t i;
    for (i = 0; i < SIZE(C); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_value_approx(ctx, COEFF(C, i), m, &tmp1);
        rational_interval_pow(&tmp2, &x_value, i);
        // tracef("tmp2 = "); lp_rational_interval_print(&tmp2, trace_out); tracef("\n");
        // tracef("tmp1 = "); lp_rational_interval_print(&tmp1, trace_out); tracef("\n");
        rational_interval_mul(&tmp2, &tmp2, &tmp1);
        // tracef("tmp2 = "); lp_rational_interval_print(&tmp2, trace_out); tracef("\n");
        // tracef("result = "); lp_rational_interval_print(&result, trace_out); tracef("\n");
        rational_interval_add(&result, &result, &tmp2);
        // tracef("result = "); lp_rational_interval_print(&result, trace_out); tracef("\n");
      }
    }

    lp_rational_interval_swap(&result, value);
    lp_rational_interval_destruct(&x_value);
    lp_rational_interval_destruct(&tmp1);
    lp_rational_interval_destruct(&tmp2);
    lp_rational_interval_destruct(&result);
  }

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_value_approx() => "); lp_rational_interval_print(value, trace_out); tracef("\n");
  }

}


/**
 * C is an univariate polynomial C(x), we compute a bound L = 1/2^k such that
 * any root of C(x) that is not zero is outsize of [L, -L].
 *
 * If C(x) = c_n x^n + ... + c_1 x + c_0 then the lower bound of on roots is an
 * upper bound 1/L on roots of C(1/x) = c_n + ... + c0 x^n.
 *
 * We compute the upper bound using the Cauchy bound
 *
 *   bound = max(|c_0|, ... |c_n|)/c_0
 *
 * We take L = 1/bound, and so we compute k = log(max/c0) = log(max) - log(c0)
 */
unsigned coefficient_root_lower_bound(const coefficient_t* C) {

  assert(C->type == COEFFICIENT_POLYNOMIAL);
  assert(coefficient_is_univariate(C));

  size_t i = 0;

  // Find the first non-zero coefficient
  while (integer_is_zero(lp_Z, &COEFF(C, i)->value.num)) {
    ++ i;
    assert(i < SIZE(C));
  }

  // First one (modulo the intial zeroes
  unsigned log_c0 = integer_log2_abs(&COEFF(C, i)->value.num);

  // Get thge max log
  unsigned max_log = log_c0;
  for (++ i; i < SIZE(C); ++ i) {
    assert(COEFF(C, i)->type == COEFFICIENT_NUMERIC);
    if (!integer_is_zero(lp_Z, &COEFF(C, i)->value.num)) {
      unsigned current_log = integer_log2_abs(&COEFF(C, i)->value.num);
      if (current_log > max_log) {
        max_log = current_log;
      }
    }
  }

  // Return the bound
  return max_log - log_c0;
}

STAT_DECLARE(int, coefficient, sgn)

int coefficient_sgn(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* m) {

  if (trace_is_enabled("coefficient::sgn")) {
    tracef("coefficient_sgn("); coefficient_print(ctx, C, trace_out); tracef("\n");
  }
  STAT(coefficient, sgn) ++;

  assert(ctx->K == lp_Z);

  int sgn;

  if (C->type == COEFFICIENT_NUMERIC) {
    // For numeric coefficients we're done
    sgn = integer_sgn(lp_Z, &C->value.num);
    if (trace_is_enabled("coefficient::sgn")) {
      tracef("coefficient_sgn(): constant => %d\n", sgn);
    }
  } else {
    assert(C->type == COEFFICIENT_POLYNOMIAL);

    if (trace_is_enabled("coefficient::sgn")) {
      tracef("coefficient_sgn(): evaluating in rationals\n");
    }

    // Try to evaluate in the rationals
    coefficient_t C_rat;
    coefficient_construct(ctx, &C_rat);
    lp_integer_t multiplier;
    integer_construct(&multiplier);
    coefficient_evaluate_rationals(ctx, C, m, &C_rat, &multiplier);

    if (trace_is_enabled("coefficient::sgn")) {
      tracef("coefficient_sgn(): C_rat = "); coefficient_print(ctx, &C_rat, trace_out); tracef("\n");
    }

    // If constant, we're done
    if (C_rat.type == COEFFICIENT_NUMERIC) {
      // val(C) = C_rat/multiplier with multiplier positive
      sgn = integer_sgn(lp_Z, &C_rat.value.num);
      if (trace_is_enabled("coefficient::sgn")) {
        tracef("coefficient_sgn(): constant => %d\n", sgn);
      }
    } else {

      // Approximate the value of C_rat
      lp_rational_interval_t C_rat_approx;
      lp_rational_interval_construct_zero(&C_rat_approx);

      if (trace_is_enabled("coefficient::sgn")) {
        tracef("coefficient_sgn(): approximating with intervals\n");
      }

      // Approximate the value by doing interval computation
      coefficient_value_approx(ctx, &C_rat, m, &C_rat_approx);

      if (trace_is_enabled("coefficient::sgn")) {
        tracef("coefficient_sgn(): approx => "); lp_rational_interval_print(&C_rat_approx, trace_out); tracef("\n");
      }

      if (C_rat_approx.is_point || !lp_rational_interval_contains_zero(&C_rat_approx)) {
        // Safe to give the sign based on the interval bound
        sgn = lp_rational_interval_sgn(&C_rat_approx);
        if (trace_is_enabled("coefficient::sgn")) {
          tracef("coefficient_sgn(): interval is good => %d\n", sgn);
        }
      } else {

        //
        // At this point the value is most likely 0.
        //

        //
        // We're still not sure, we need to evaluate the sign using resultants.
        // We construct the polynomial
        //
        //   A(z, x1, ..., xn) = z - C_rat(x1, ..., xn).
        //
        // If the value of C_rat is a, and values of x_i is a_i with defining polynomial p_i(x_i)
        // then we know that
        //
        //   B(z) = Resultant(A, p1, ..., pn)
        //
        // has a as at least one zero.
        //
        // We can estimate the bound L on the smallest root of B, and refine the
        // evaluation interval until it is either in (-L, L), in which case
        // it must be 0, or until it doesn't contain 0, in which case we get
        // the sign.
        //

        // The temporary variable, we'll be using
        lp_variable_t z = lp_polynomial_context_get_temp_variable(ctx);

        // A = z - C_rat
        coefficient_t A;
        lp_integer_t A_lc;
        integer_construct_from_int(lp_Z, &A_lc, 1);
        coefficient_construct_simple(ctx, &A, &A_lc, z, 1);
        coefficient_sub(ctx, &A, &A, &C_rat);
        lp_integer_destruct(&A_lc);

        // List of variables in C_rat;
        lp_variable_list_t C_rat_vars;
        lp_variable_list_construct(&C_rat_vars);
        coefficient_get_variables(&C_rat, &C_rat_vars);

        // Cache the assignment intervals
        lp_dyadic_interval_t* rational_interval_cache = algebraic_interval_remember(&C_rat_vars, m);

        // Compute B (we keep it in A) by resolving out all the variables except z
        lp_variable_order_make_bot(ctx->var_order, z);
        coefficient_order(ctx, &A);
        coefficient_resolve_algebraic(ctx, &A, m, &A);
        lp_variable_order_make_bot(ctx->var_order, lp_variable_null);

        // Get the lower bound on the roots
        unsigned k = coefficient_root_lower_bound(&A);
        // Interval (-L, L)
        lp_rational_interval_t L_interval;
        lp_rational_t L, L_neg;
        rational_construct_from_int(&L, 1, 1);
        rational_div_2exp(&L, &L, k);
        rational_construct(&L_neg);
        rational_neg(&L_neg, &L);
        lp_rational_interval_construct(&L_interval, &L_neg, 1, &L, 1);

        // Refine until done, i.e. C_rat_approx = (l, u) not contains 0, or
        //
        for (;;) {

          if (trace_is_enabled("coefficient::sgn")) {
            tracef("coefficient_sgn(): C_rat_approx = "); lp_rational_interval_print(&C_rat_approx, trace_out); tracef("\n");
            tracef("coefficient_sgn(): L_interval = "); lp_rational_interval_print(&L_interval, trace_out); tracef("\n");
          }

          // If a point, we're done
          if (lp_rational_interval_is_point(&C_rat_approx)) {
            break;
          }

          // If no zero in interval, we're done
          if (!lp_rational_interval_contains_zero(&C_rat_approx)) {
            break;
          }

          // If contained in L, we're also done
          int contains_a = lp_rational_interval_contains_rational(&L_interval, &C_rat_approx.a);
          int contains_b = lp_rational_interval_contains_rational(&L_interval, &C_rat_approx.b);
          if (contains_a && contains_b) {
            assert(lp_rational_interval_contains_zero(&C_rat_approx));
            break;
          }

          // Refine the values
          size_t i;
          for (i = 0; i < C_rat_vars.list_size; ++ i) {
            lp_variable_t x_i = C_rat_vars.list[i];
            const lp_value_t* x_i_value = lp_assignment_get_value(m, x_i);
            if (!lp_value_is_rational(x_i_value)) {
              lp_algebraic_number_refine_const(&x_i_value->value.a);
            }
          }

          // Approximate the value by doing interval computation
          coefficient_value_approx(ctx, C, m, &C_rat_approx);
        }

        // Restore the cache (will free the cache)
        algebraic_interval_restore(&C_rat_vars, rational_interval_cache, m);

        // Release the temporary variable
        lp_polynomial_context_release_temp_variable(ctx, z);

        // Destruct temps
        coefficient_destruct(&A);
        lp_variable_list_destruct(&C_rat_vars);
        lp_rational_interval_destruct(&L_interval);
        lp_rational_destruct(&L);
        lp_rational_destruct(&L_neg);

        // Safe to give the sign based on the interval bound
        sgn = lp_rational_interval_sgn(&C_rat_approx);
        if (trace_is_enabled("coefficient::sgn")) {
          tracef("coefficient_sgn(): interval is good => %d\n", sgn);
        }
      }

      // Destruct temps
      lp_rational_interval_destruct(&C_rat_approx);
    }

    // Destruct temps
    lp_integer_destruct(&multiplier);
    coefficient_destruct(&C_rat);
  }

  if (trace_is_enabled("coefficient::sgn")) {
    tracef("coefficient_sgn() => %d\n", sgn);
  }

  if (sgn < 0) { return -1; }
  if (sgn > 0) { return  1; }

  return 0;
}

STAT_DECLARE(int, coefficient, lc_sgn)

int coefficient_lc_sgn(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  STAT(coefficient, lc_sgn) ++;

  while (C->type != COEFFICIENT_NUMERIC) {
    C = coefficient_lc(C);
  }

  return integer_sgn(ctx->K, &C->value.num);
}


STAT_DECLARE(int, coefficient, in_order)

int coefficient_in_order(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  TRACE("coefficient::internal", "coefficient_in_order()\n");
  STAT(coefficient, in_order) ++;

  size_t i;
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
        if (lp_variable_order_cmp(ctx->var_order, VAR(C), VAR(C_i)) <= 0) {
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


int coefficient_cmp_general(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2, int compare_values) {
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
    int var_cmp = lp_variable_order_cmp(ctx->var_order, VAR(C1), VAR(C2));
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

STAT_DECLARE(int, coefficient, cmp)

int coefficient_cmp(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_cmp()\n");
  STAT(coefficient, cmp) ++;

  return coefficient_cmp_general(ctx, C1, C2, 1);
}

STAT_DECLARE(int, coefficient, cmp_type)

int coefficient_cmp_type(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::internal", "coefficient_cmp_type()\n");
  STAT(coefficient, cmp_type) ++;

  return coefficient_cmp_general(ctx, C1, C2, 0);
}

void coefficient_traverse(const lp_polynomial_context_t* ctx, const coefficient_t* C, traverse_f f, lp_monomial_t* m, void* data) {

  if (trace_is_enabled("coefficient::order")) {
    tracef("order = "); lp_variable_order_print(ctx->var_order, ctx->var_db, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
    tracef("m = "); monomial_print(ctx, m, trace_out); tracef("\n");
  }

  size_t d;
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
        lp_monomial_push(m, VAR(C), d);
        coefficient_traverse(ctx, COEFF(C, d), f, m, data);
        lp_monomial_pop(m);
      }
    }
    break;
  }
}

/**
 * Method called to add a monomial to C. The monomial should be ordered in the
 * same order as C, top variable at the m[0].
 */
void coefficient_add_ordered_monomial(const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* C_void) {

  coefficient_t* C = (coefficient_t*) C_void;

  if (trace_is_enabled("coefficient::order")) {
    tracef("coefficient_add_monomial():\n");
    tracef("m = "); monomial_print(ctx, m, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
  }

  if (m->n == 0) {
    // Just add a constant to C
    switch (C->type) {
    case COEFFICIENT_NUMERIC:
      integer_add(ctx->K, &C->value.num, &C->value.num, &m->a);
      break;
    case COEFFICIENT_POLYNOMIAL:
      coefficient_add_ordered_monomial(ctx, m, COEFF(C, 0));
      break;
    }
  } else {
    // Proper monomial (first variable is the top one)
    lp_variable_t x = m->p[0].x;
    unsigned d = m->p[0].d;
    // Compare the variables
    if (C->type == COEFFICIENT_NUMERIC || lp_variable_order_cmp(ctx->var_order, x, VAR(C)) >= 0) {
      coefficient_ensure_capacity(ctx, C, x, d+1);
      // Now, add the monomial to the right place
      m->p ++;
      m->n --;
      coefficient_add_ordered_monomial(ctx, m, COEFF(C, d));
      coefficient_normalize(ctx, C);
      m->p --;
      m->n ++;
    } else {
      coefficient_add_ordered_monomial(ctx, m, COEFF(C, 0));
    }
  }

  assert(coefficient_is_normalized(ctx, C));
}

void coefficient_order_and_add_monomial(const lp_polynomial_context_t* ctx, lp_monomial_t* m, void* C_void) {
  lp_monomial_t m_ordered;
  lp_monomial_construct_copy(ctx, &m_ordered, m, /** sort */ 1);
  coefficient_add_ordered_monomial(ctx, &m_ordered, C_void);
  lp_monomial_destruct(&m_ordered);
}

void coefficient_add_monomial(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_monomial_t* m) {
  lp_monomial_t m_ordered;
  lp_monomial_construct_copy(ctx, &m_ordered, m, /** sort */ 1);
  coefficient_add_ordered_monomial(ctx, &m_ordered, C);
  lp_monomial_destruct(&m_ordered);
}

STAT_DECLARE(int, coefficient, order)

void coefficient_order(const lp_polynomial_context_t* ctx, coefficient_t* C) {
  TRACE("coefficient", "coefficient_order()\n");
  STAT(coefficient, order) ++;

  if (C->type == COEFFICIENT_NUMERIC) {
    // Numeric coefficients are always OK
    return;
  }

  if (trace_is_enabled("coefficient::order")) {
    tracef("order = "); lp_variable_order_print(ctx->var_order, ctx->var_db, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
  }

  // The coefficient we are building
  coefficient_t result;
  coefficient_construct(ctx, &result);
  // The monomials build in the original order
  lp_monomial_t m_tmp;
  lp_monomial_construct(ctx, &m_tmp);
  // For each monomial of C, add it to the result
  coefficient_traverse(ctx, C, coefficient_order_and_add_monomial, &m_tmp, &result);
  // Keep the result
  coefficient_swap(C, &result);
  // Destroy temps
  lp_monomial_destruct(&m_tmp);
  coefficient_destruct(&result);

  assert(coefficient_is_normalized(ctx, C));
}

STAT_DECLARE(int, coefficient, add)

#define MAX(x, y) (x >= y ? x : y)

void coefficient_add(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_add()\n");
  STAT(coefficient, add) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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
      size_t i;
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

  if (trace_is_enabled("coefficient::arith")) {
    tracef("add = "); coefficient_print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, neg)

void coefficient_neg(const lp_polynomial_context_t* ctx, coefficient_t* N, const coefficient_t* C) {
  TRACE("coefficient::arith", "coefficient_neg()\n");
  STAT(coefficient, neg) ++;

  size_t i;
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

STAT_DECLARE(int, coefficient, sub)

void coefficient_sub(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_sub()\n");
  STAT(coefficient, sub) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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
      size_t i;
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


void coefficient_add_mul(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2);

STAT_DECLARE(int, coefficient, mul)

void coefficient_mul(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_mul()\n");
  STAT(coefficient, mul) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
  }

  size_t i, j;
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
              if (trace_is_enabled("coefficient::arith")) {
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
    // C1 < C2, multiply each coefficient of C2 with C1
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

  if (trace_is_enabled("coefficient::arith")) {
    tracef("mul = "); coefficient_print(ctx, P, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, P));
}

STAT_DECLARE(int, coefficient, mul_int)

void coefficient_mul_int(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, long a) {
  TRACE("coefficient::arith", "coefficient_mul_int()\n");
  STAT(coefficient, mul_int) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
    tracef("n  = %ld\n", a);
  }

  size_t i;
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

  if (trace_is_enabled("coefficient::arith")) {
    tracef("mul = "); coefficient_print(ctx, P, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, P));
}

STAT_DECLARE(int, coefficient, mul_integer)

void coefficient_mul_integer(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, const lp_integer_t* a) {
  TRACE("coefficient::arith", "coefficient_mul_int()\n");
  STAT(coefficient, mul_int) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
    tracef("n  = "); integer_print(a, trace_out);
  }

  size_t i;
  coefficient_t result;

  if (C->type == COEFFICIENT_NUMERIC) {
    if (P->type == COEFFICIENT_POLYNOMIAL) {
      coefficient_construct(ctx, &result);
      integer_mul(ctx->K, &result.value.num, &C->value.num, a);
      coefficient_swap(&result, P);
      coefficient_destruct(&result);
    } else {
      integer_mul(ctx->K, &P->value.num, &C->value.num, a);
    }
  } else {
    coefficient_construct_rec(ctx, &result, VAR(C), SIZE(C));
    for (i = 0; i < SIZE(C); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        coefficient_mul_integer(ctx, COEFF(&result, i), COEFF(C, i), a);
      }
    }
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, P);
    coefficient_destruct(&result);
  }

  assert(coefficient_is_normalized(ctx, P));
}


STAT_DECLARE(int, coefficient, shl)

void coefficient_shl(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, lp_variable_t x, unsigned n) {
  TRACE("coefficient::arith", "coefficient_shl()\n");
  STAT(coefficient, shl) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
    tracef("x = %s\n", lp_variable_db_get_name(ctx->var_db, x));
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

  if (trace_is_enabled("coefficient::arith")) {
    tracef("coefficient_shl() =>"); coefficient_print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, shr)

void coefficient_shr(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C, unsigned n) {
  TRACE("coefficient::arith", "coefficient_shr()\n");
  STAT(coefficient, shl) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
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
    for (i = 0; i < (int) SIZE(C) - (int) n; ++ i) {
      coefficient_assign(ctx, COEFF(&result, i), COEFF(C, i + n));
    }
    coefficient_swap(&result, S);
    coefficient_destruct(&result);
  }

  if (trace_is_enabled("coefficient::arith")) {
    tracef("coefficient_shr() =>"); coefficient_print(ctx, S, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, S));
}

STAT_DECLARE(int, coefficient, pow)

void coefficient_pow(const lp_polynomial_context_t* ctx, coefficient_t* P, const coefficient_t* C, unsigned n) {
  TRACE("coefficient", "coefficient_pow()\n");
  STAT(coefficient, pow) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("P = "); coefficient_print(ctx, P, trace_out); tracef("\n");
    tracef("C = "); coefficient_print(ctx, C, trace_out); tracef("\n");
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

STAT_DECLARE(int, coefficient, add_mul)

void coefficient_add_mul(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient::arith", "coefficient_add_mul()\n");
  STAT(coefficient, add_mul) ++;

  if (trace_is_enabled("coefficient::arith")) {
    tracef("S = "); coefficient_print(ctx, S, trace_out); tracef("\n");
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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

STAT_DECLARE(int, coefficient, sub_mul)

void coefficient_sub_mul(const lp_polynomial_context_t* ctx, coefficient_t* S, const coefficient_t* C1, const coefficient_t* C2) {
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

STAT_DECLARE(int, coefficient, derivative)

void coefficient_derivative(const lp_polynomial_context_t* ctx, coefficient_t* C_d, const coefficient_t* C) {
  TRACE("coefficient", "coefficient_derivative()\n");
  STAT(coefficient, derivative) ++;

  size_t i;
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

void coefficient_div_degrees(const lp_polynomial_context_t* ctx, coefficient_t* C, size_t p) {
  if (C->type == COEFFICIENT_POLYNOMIAL) {
    size_t i;
    for (i = 1; i < SIZE(C); ++ i) {
      if (!coefficient_is_zero(ctx, COEFF(C, i))) {
        assert(i % p == 0);
        assert(coefficient_is_zero(ctx, COEFF(C, i/ p)));
        coefficient_swap(COEFF(C, i), COEFF(C, i / p));
      }
    }
    coefficient_normalize(ctx, C);
  }
}

///
/// Forward declarations of division/reduction/gcd stuff
///

void coefficient_div(const lp_polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2);

//
// Implementation of the division/reduction/gcd stuff
//

STAT_DECLARE(int, coefficient, reduce)

void coefficient_reduce(
    const lp_polynomial_context_t* ctx,
    const coefficient_t* A, const coefficient_t* B,
    coefficient_t* P, coefficient_t* Q, coefficient_t* R,
    remaindering_type_t type)
{
  TRACE("coefficient", "coefficient_reduce()\n");
  STAT(coefficient, reduce) ++;

  if (trace_is_enabled("coefficient::reduce")) {
    tracef("A = "); coefficient_print(ctx, A, trace_out); tracef("\n");
    tracef("B = "); coefficient_print(ctx, B, trace_out); tracef("\n");
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
  lp_variable_t x = VAR(A);

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

    if (trace_is_enabled("coefficient::reduce")) {
      tracef("lcm = "); coefficient_print(ctx, &lcm, trace_out); tracef("\n");
      tracef("R = "); coefficient_print(ctx, &R_tmp, trace_out); tracef("\n");
      tracef("r = "); coefficient_print(ctx, &r, trace_out); tracef("\n");
      tracef("B = "); coefficient_print(ctx, B, trace_out); tracef("\n");
      tracef("b = "); coefficient_print(ctx, &b, trace_out); tracef("\n");
    }

    // R'
    coefficient_mul(ctx, &R_tmp, &R_tmp, &r);
    coefficient_sub_mul(ctx, &R_tmp, B, &b);

    if (trace_is_enabled("coefficient::reduce")) {
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

    if (trace_is_enabled("coefficient::reduce")) {
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


STAT_DECLARE(int, coefficient, div)

void coefficient_div(const lp_polynomial_context_t* ctx, coefficient_t* D, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_div()\n");
  STAT(coefficient, div) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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
    if (trace_is_enabled("coefficient::check_division")) {
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_div() => "); coefficient_print(ctx, D, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, D));
}

STAT_DECLARE(int, coefficient, rem)

void coefficient_rem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_rem()\n");
  STAT(coefficient, rem) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_rem() => "); coefficient_print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, sprem)

void coefficient_sprem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_sprem()\n");
  STAT(coefficient, sprem) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_sprem() => "); coefficient_print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, prem)

void coefficient_prem(const lp_polynomial_context_t* ctx, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_prem()\n");
  STAT(coefficient, rem) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_prem() => "); coefficient_print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}


STAT_DECLARE(int, coefficient, divrem)

void coefficient_divrem(const lp_polynomial_context_t* ctx, coefficient_t* D, coefficient_t* R, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_divrem()\n");
  STAT(coefficient, divrem) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("C1 = "); coefficient_print(ctx, C1, trace_out); tracef("\n");
    tracef("C2 = "); coefficient_print(ctx, C2, trace_out); tracef("\n");
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_divrem() => \n");
    tracef("D = "); coefficient_print(ctx, D, trace_out); tracef("\n");
    tracef("D = "); coefficient_print(ctx, R, trace_out); tracef("\n");
  }

  assert(coefficient_is_normalized(ctx, R));
}

STAT_DECLARE(int, coefficient, divides)

int coefficient_divides(const lp_polynomial_context_t* ctx, const coefficient_t* C1, const coefficient_t* C2) {
  TRACE("coefficient", "coefficient_divides()\n");
  STAT(coefficient, divides) ++;

  coefficient_t R;
  coefficient_construct(ctx, &R);
  coefficient_prem(ctx, &R, C2, C1);
  int divides = coefficient_is_zero(ctx, &R);
  coefficient_destruct(&R);

  return divides;
}

int coefficient_is_univariate(const coefficient_t* C) {
  size_t i;
  if (C->type == COEFFICIENT_NUMERIC) {
    return 1;
  } else {
    for (i = 0; i < SIZE(C); ++ i) {
      if (COEFF(C, i)->type != COEFFICIENT_NUMERIC) {
        return 0;
      }
    }
    return 1;
  }
}

const lp_integer_t* coefficient_get_constant(const coefficient_t* C) {
  while (C->type == COEFFICIENT_POLYNOMIAL) {
    C = COEFF(C, 0);
  }
  return &C->value.num;
}

lp_upolynomial_t* coefficient_to_univariate(const lp_polynomial_context_t* ctx, const coefficient_t* C) {
  assert(C->type == COEFFICIENT_POLYNOMIAL);

  lp_integer_t* coeff = malloc(sizeof(lp_integer_t)*SIZE(C));

  size_t i;
  for (i = 0; i < SIZE(C); ++ i) {
    integer_construct_copy(ctx->K, coeff + i, coefficient_get_constant(COEFF(C, i)));
  }

  lp_upolynomial_t* C_u = lp_upolynomial_construct(ctx->K, SIZE(C) - 1, coeff);

  for (i = 0; i < SIZE(C); ++ i) {
    integer_destruct(coeff + i);
  }
  free(coeff);

  return C_u;
}

STAT_DECLARE(int, coefficient, resultant)

void coefficient_resultant(const lp_polynomial_context_t* ctx, coefficient_t* res, const coefficient_t* A, const coefficient_t* B) {

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_resultant("); coefficient_print(ctx, A, trace_out); tracef(", "); coefficient_print(ctx, B, trace_out); tracef(")\n");
  }
  STAT(coefficient, resultant) ++;

  if (trace_is_enabled("coefficient")) {
    tracef("A = "); coefficient_print(ctx, A, trace_out); tracef("\n");
    tracef("B = "); coefficient_print(ctx, B, trace_out); tracef("\n");
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
  size_t i;
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

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_resultant() => "); coefficient_print(ctx, res, trace_out); tracef("\n");
  }

}

///
/// Normalization that everyone is using
///

static void
coefficient_normalize(const lp_polynomial_context_t* ctx, coefficient_t* C) {
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

static void
coefficient_normalize_m(const lp_polynomial_context_t* ctx, coefficient_t* C, const lp_assignment_t* m) {
  if (C->type == COEFFICIENT_POLYNOMIAL) {
    assert(C->value.rec.size >= 1);
    size_t i = C->value.rec.size - 1;
    // Find the first non-zero coefficient
    while (i > 0 && coefficient_sgn(ctx, COEFF(C, i), m) == 0) {
      i --;
      coefficient_assign_int(ctx, COEFF(C, i), 0);
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
coefficient_is_normalized(const lp_polynomial_context_t* ctx, coefficient_t* C) {
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
coefficient_ensure_capacity(const lp_polynomial_context_t* ctx, coefficient_t* C, lp_variable_t x, size_t capacity) {
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
      assert(lp_variable_order_cmp(ctx->var_order, x, VAR(C)) > 0);
      // Same as for constants above
      coefficient_construct_rec(ctx, &tmp, x, capacity);
      coefficient_swap(COEFF(&tmp, 0), C);
      coefficient_swap(C, &tmp);
      coefficient_destruct(&tmp);
    } else if (capacity > C->value.rec.capacity) {
      // Already recursive polynomial, resize up
      C->value.rec.coefficients = realloc(C->value.rec.coefficients, capacity * sizeof(coefficient_t));
      size_t i;
      for (i = C->value.rec.capacity; i < capacity; ++ i) {
        coefficient_construct(ctx, C->value.rec.coefficients + i);
      }
      C->value.rec.capacity  = capacity;
      C->value.rec.size = capacity;
    }
    break;
  }
}

void coefficient_evaluate_rationals(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M, coefficient_t* C_out, lp_integer_t* multiplier) {

  assert(multiplier);
  assert(ctx->K == lp_Z);

  size_t i;
  lp_variable_t x;

  // Start wit multiplier 1
  integer_assign_int(lp_Z, multiplier, 1);

  if (C->type == COEFFICIENT_NUMERIC) {
    // Just a number, we're done
    coefficient_assign(ctx, C_out, C);
  } else {
    assert(C->type == COEFFICIENT_POLYNOMIAL);

    // Temp for the result
    coefficient_t result;

    // Get the variable and it's value, if any
    x = VAR(C);
    const lp_value_t* x_value = lp_assignment_get_value(M, x);

    // The degree of the polynomial
    size_t size = SIZE(C);

    // Check if the value is rational and we can substitute it
    if (!lp_value_is_rational(x_value))
    {
      // We can not substitute so
      //
      //   C = a_n * x^n + ... + a_1 * x + a_0
      //
      // We substitute in all a_n obtaining a_k = b_n / m_k, m = lcm(m_1, ..., m_n)
      //
      //   m * c = sum     b_k * x^k * m / m_k

      coefficient_construct_rec(ctx, &result, x, size);

      // Compute the evaluation of the coefficients
      lp_integer_t* m = malloc(sizeof(coefficient_t)*size);
      for (i = 0; i < size; ++ i) {
        integer_construct(m + i);
        coefficient_evaluate_rationals(ctx, COEFF(C, i), M, COEFF(&result, i), m + i);
      }

      // Compute the lcm of the m's
      lp_integer_construct_copy(lp_Z, multiplier, m);
      for (i = 1; i < size; ++ i) {
        integer_lcm_Z(multiplier, multiplier, m + i);
      }

      // Sum up
      lp_integer_t tmp;
      integer_construct(&tmp);
      for (i = 0; i < size; ++ i) {
        // m / m_k
        integer_div_exact(lp_Z, &tmp, multiplier, m + i);
        // b_i = b_i * R
        coefficient_mul_integer(ctx, COEFF(&result, i), COEFF(&result, i), &tmp);
      }
      integer_destruct(&tmp);

      // Remove the temps
      for (i = 0; i < size; ++ i) {
        integer_destruct(m + i);
      }
      free(m);

    } else {

      coefficient_construct(ctx, &result);

      // We have a value value = p/q
      lp_integer_t p, q;
      integer_construct(&p);
      integer_construct(&q);
      lp_value_get_num(x_value, &p);
      lp_value_get_den(x_value, &q);

      // If we can substitute then
      //
      //   C = a_n * (p/q)^n + ... + a_1 * (p/q) + a_0
      //
      // We substitute in all a_n obtaining a_k = b_n / m_k and get
      //
      //   q^n * C = b_n * (p^n/m_n) + ... + b_1 * (p*q^n-1/m_1) + b_0 * q^n/m_0
      //
      // We get the m = lcm(m_1, ..., m_n) and get
      //
      //   q^n * m * c = sum     b_k * p^k * q^(n-k) * m / m_k


      // Compute the evaluation of the coefficients
      coefficient_t* b = malloc(sizeof(coefficient_t)*size);
      lp_integer_t* m = malloc(sizeof(lp_integer_t)*size);
      for (i = 0; i < size; ++ i) {
        coefficient_construct(ctx, b + i);
        integer_construct(m + i);
        coefficient_evaluate_rationals(ctx, COEFF(C, i), M, b + i, m + i);
      }

      // Compute the lcm of the m's
      lp_integer_t m_lcm;
      lp_integer_construct_copy(lp_Z, &m_lcm, m);
      for (i = 1; i < size; ++ i) {
        integer_lcm_Z(&m_lcm, &m_lcm, m + i);
      }

      // The powers
      lp_integer_t p_power, q_power;
      integer_construct_from_int(lp_Z, &p_power, 1);
      integer_construct(&q_power);
      integer_pow(lp_Z, &q_power, &q, size-1);

      // Set the multiplier
      integer_mul(lp_Z, multiplier, &q_power, &m_lcm);

      // Sum up
      lp_integer_t R;
      integer_construct(&R);
      for (i = 0; i < size; ++ i) {
        if (i) {
          // Update powers
          integer_mul(lp_Z, &p_power, &p_power, &p);
          integer_div_exact(lp_Z, &q_power, &q_power, &q);
        }
        // R = p^i * q^(n-i) * m / m_k
        integer_div_exact(lp_Z, &R, &m_lcm, m + i);
        integer_mul(lp_Z, &R, &R, &p_power);
        integer_mul(lp_Z, &R, &R, &q_power);
        // b_i = b_i * R
        coefficient_mul_integer(ctx, b + i, b + i, &R);
        // Add it
        coefficient_add(ctx, &result, &result, b + i);
      }
      integer_destruct(&R);

      // Remove the temps
      for (i = 0; i < size; ++ i) {
        coefficient_destruct(b + i);
        integer_destruct(m + i);
      }
      free(b);
      free(m);
      integer_destruct(&m_lcm);
      integer_destruct(&p);
      integer_destruct(&q);
      integer_destruct(&p_power);
      integer_destruct(&q_power);
    }

    // Finish up
    coefficient_normalize(ctx, &result);
    coefficient_swap(&result, C_out);
    coefficient_destruct(&result);
  }

  assert(integer_sgn(lp_Z, multiplier) > 0);

}

void coefficient_get_variables(const coefficient_t* C, lp_variable_list_t* vars) {
  if (C->type != COEFFICIENT_NUMERIC) {
    // Add the variable, if not already there
    lp_variable_t x = VAR(C);
    if (lp_variable_list_index(vars, x) < 0) {
      lp_variable_list_push(vars, x);
    }
    // Add children
    size_t i;
    for (i = 0; i < SIZE(C); ++ i) {
      coefficient_get_variables(COEFF(C, i), vars);
    }
  }
}

/**
 * Isolate out the roots of a univariate polynomial.
 */
static
void coefficient_roots_isolate_univariate(const lp_polynomial_context_t* ctx, const coefficient_t* A, lp_value_t* roots, size_t* roots_size) {

  if (trace_is_enabled("coefficient::roots")) {
    tracef("coefficient_roots_isolate(): univariate, root finding\n");
    tracef("coefficient_roots_isolate(): A = "); coefficient_print(ctx, A, trace_out); tracef("\n");
  }

  // Special case for linear polynomials
  if (coefficient_degree(A) == 1) {
    // A_rat = ax + b => root = -b/a
    *roots_size = 1;
    lp_rational_t root;
    rational_construct_from_div(&root, &COEFF(A, 0)->value.num, &COEFF(A, 1)->value.num);
    rational_neg(&root, &root);
    lp_value_construct(roots, LP_VALUE_RATIONAL, &root);
    rational_destruct(&root);
  } else {
    // Get the univariate version
    lp_upolynomial_t* A_u = coefficient_to_univariate(ctx, A);
    // Make space for the algebraic numbers
    lp_algebraic_number_t* algebraic_roots = malloc(lp_upolynomial_degree(A_u) * sizeof(lp_algebraic_number_t));

    if (trace_is_enabled("coefficient::roots")) {
      tracef("coefficient_roots_isolate(): A_u = "); lp_upolynomial_print(A_u, trace_out); tracef("\n");
    }
    // Isolate
    lp_upolynomial_roots_isolate(A_u, algebraic_roots, roots_size);

    // Copy over the roots
    size_t i;
    for (i = 0; i < *roots_size; ++ i) {
      lp_value_construct(roots + i, LP_VALUE_ALGEBRAIC, algebraic_roots + i);
      lp_algebraic_number_destruct(algebraic_roots + i);
    }
    // Free the temps
    free(algebraic_roots);
    lp_upolynomial_delete(A_u);
  }
}

/**
 * Resolve out the algebraic numbers assignment.
 */
static
void coefficient_resolve_algebraic(const lp_polynomial_context_t* ctx, const coefficient_t* A, const lp_assignment_t* m, coefficient_t* A_alg) {

  coefficient_assign(ctx, A_alg, A);

  // Eliminate all variables
  for (;;) {

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_resolve_algebraic(): A_alg = "); coefficient_print(ctx, A_alg, trace_out); tracef("\n");
    }

    // If no more variables, it either has no zeroes, or it vanished
    if (A_alg->type == COEFFICIENT_NUMERIC) {
      break;
    }

    // The variable to eliminate
    lp_variable_t y = VAR(A_alg);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_resolve_algebraic(): resolving out %s\n", lp_variable_db_get_name(ctx->var_db, y));
    }

    // Get the algebraic number of the value of y as a polynomial in y
    const lp_value_t* y_value = lp_assignment_get_value(m, y);
    if (y_value->type == LP_VALUE_NONE) {
      assert(coefficient_is_univariate(A_alg));
      break;
    }

    assert(!lp_value_is_rational(y_value));
    const lp_upolynomial_t* y_upoly = y_value->value.a.f;
    assert(y_upoly != 0);

    // Make it multivariate
    coefficient_t y_poly;
    coefficient_construct_from_univariate(ctx, &y_poly, y_upoly, y);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_resolve_algebraic(): y_poly = "); coefficient_print(ctx, &y_poly, trace_out); tracef("\n");
    }

    // Resolve
    coefficient_resultant(ctx, A_alg, A_alg, &y_poly);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_resolve_algebraic(): resultant done = "); coefficient_print(ctx, A_alg, trace_out); tracef("\n");
    }

    // Destroy temp
    coefficient_destruct(&y_poly);
  }
}


void coefficient_roots_isolate(const lp_polynomial_context_t* ctx, const coefficient_t* A, const lp_assignment_t* M, lp_value_t* roots, size_t* roots_size) {

  if (trace_is_enabled("coefficient::roots")) {
    tracef("coefficient_roots_isolate("); coefficient_print(ctx, A, trace_out); tracef(")\n");
  }

  assert(A->type == COEFFICIENT_POLYNOMIAL);
  assert(lp_assignment_get_value(M, coefficient_top_variable(A))->type == LP_VALUE_NONE);

  // Evaluate in the rationals
  coefficient_t A_rat;
  lp_integer_t multiplier;
  integer_construct(&multiplier);
  coefficient_construct(ctx, &A_rat);
  coefficient_evaluate_rationals(ctx, A, M, &A_rat, &multiplier);

  if (trace_is_enabled("coefficient::roots")) {
    tracef("coefficient_roots_isolate(): rational evaluation: "); coefficient_print(ctx, &A_rat, trace_out); tracef("\n");
  }

  // If this is a constant polynomial, no zeroes
  if (coefficient_is_constant(&A_rat)) {
    if (trace_is_enabled("coefficient::roots")) {
      tracef("coefficient_roots_isolate(): constant, no roots\n");
    }
    *roots_size = 0;
  } else {

    // TODO: normalize A_rat by gcd of coefficients

    // If univariate, just isolate the univariate roots
    if (coefficient_is_univariate(&A_rat)) {
      coefficient_roots_isolate_univariate(ctx, &A_rat, roots, roots_size);
    } else {

      // Can not evaluate in Q, do the expensive stuff
      if (trace_is_enabled("coefficient::roots")) {
        tracef("coefficient_roots_isolate(): not univariate\n");
      }

      // Not univariate, A_rat(x, y0, ..., y_n), with x the top variable. We
      // eliminate the y variables with resolvents of the algebraic numbers to
      // obtain a univariate polynomial that includes the zeroes of A_rat.
      // We reverse the order, and just eliminate one-by one, until we reach x

      // The variable
      lp_variable_t x = VAR(&A_rat);
      assert(x == VAR(A));

      // Top variable must be unassigned
      assert(lp_assignment_get_value(M, x)->type == LP_VALUE_NONE);

      // Rerder to (y_n, ..., y0, x) and then restore the order
      lp_variable_order_make_bot(ctx->var_order, x);
      coefficient_order(ctx, &A_rat);

      // The copy where we compute the resultants
      coefficient_t A_alg;
      coefficient_construct(ctx, &A_alg);

      // Resolve algebraic assignments
      coefficient_resolve_algebraic(ctx, &A_rat, M, &A_alg);

      if (trace_is_enabled("coefficient::roots")) {
        tracef("coefficient_roots_isolate(): A_alg = "); coefficient_print(ctx, &A_alg, trace_out); tracef("\n");
      }

      if (coefficient_is_zero(ctx, &A_alg)) {

        if (trace_is_enabled("coefficient::roots")) {
          tracef("coefficient_roots_isolate(): vanished :(\n");
        }

        //
        // A_alg vanished, but this might be due to other roots of the polynomials
        // that define the assignment. We will first manually check if the polynomial
        // A_rat vanishes by trying to evaluate coefficient.
        //
        // If A_rat vanishes, we're done, there is no roots.
        //
        // If A_rat doesn't vanish, we have that (modulo the assignment)
        //
        //    A_rat = c_k x^k + ... + c_0
        //
        // If c_k evaluates to a_k, we can extend the assignment with a new
        // variable y -> a_k and get the roots of the polynomial
        //
        //   B = y x^k + ...
        //
        // This polynomial will not vanish by construction since 0 is not a zero
        // of the polynomial defining a_k (by properties of root isolation,
        // defining polynomials are don't have 0 as a root).
        //

        // Reduce based on the model
        coefficient_normalize_m(ctx, &A_rat, M);

        // If A_rat vanishes, there is not zeroes
        if (coefficient_is_zero(ctx, &A_rat)) {
          if (trace_is_enabled("coefficient::roots")) {
            tracef("coefficient_roots_isolate(): truly vanished :)\n");
          }
        } else {
          // Get the value of the coefficient
          lp_value_t* lc_value = coefficient_evaluate(ctx, coefficient_lc(&A_rat), M);

          // The new variable
          lp_variable_t y = lp_polynomial_context_get_temp_variable(ctx);

          // Set the value
          assert(lp_assignment_get_value(M, y)->type == LP_VALUE_NONE);
          lp_assignment_set_value((lp_assignment_t*) M, y, lc_value);

          // Now, do it recursively (restore the order first)
          lp_variable_order_make_bot(ctx->var_order, lp_variable_null);
          coefficient_order(ctx, &A_rat);
          coefficient_roots_isolate(ctx, &A_rat, M, roots, roots_size);

          // Undo local stuff
          lp_assignment_set_value((lp_assignment_t*) M, y, 0);
          lp_polynomial_context_release_temp_variable(ctx, y);
          lp_value_delete(lc_value);
        }
      } else if (coefficient_is_constant(&A_alg)) {
        // If constant, we have no roots
        if (trace_is_enabled("coefficient::roots")) {
          tracef("coefficient_roots_isolate(): no roots, constant\n");
        }
        *roots_size = 0;
      } else {

        // Isolate as univariate, then filter any extra roots obtained from the defining polynomials
        if (trace_is_enabled("coefficient::roots")) {
          tracef("coefficient_roots_isolate(): univariate reduct\n");
        }
        assert(coefficient_degree(&A_alg) > 1);

        // Get the univariate version
        lp_upolynomial_t* A_alg_u = coefficient_to_univariate(ctx, &A_alg);
        // Make space for the algebraic numbers
        lp_algebraic_number_t* algebraic_roots = malloc(lp_upolynomial_degree(A_alg_u) * sizeof(lp_algebraic_number_t));
        lp_upolynomial_roots_isolate(A_alg_u, algebraic_roots, roots_size);

        // We'll assigned x, so remove from bot
        lp_variable_order_make_bot(ctx->var_order, lp_variable_null);
        coefficient_order(ctx, &A_rat);

        if (trace_is_enabled("coefficient::roots")) {
          tracef("coefficient_roots_isolate(): filtering roots\n");
        }

        // Filter any bad roots
        size_t i, to_keep;
        for (i = 0, to_keep = 0; i < *roots_size; ++ i) {
          // Set the value of the variable
          assert(lp_assignment_get_value(M, x)->type == LP_VALUE_NONE);
          lp_value_t x_value;
          lp_value_construct(&x_value, LP_VALUE_ALGEBRAIC, algebraic_roots + i);
          lp_assignment_set_value((lp_assignment_t*) M, x, &x_value);

          if (trace_is_enabled("coefficient::roots")) {
            tracef("coefficient_roots_isolate(): checking root: "); lp_value_print(&x_value, trace_out); tracef("\n");
          }

          if (coefficient_sgn(ctx, &A_rat, M) == 0) {
            if (i != to_keep) {
              lp_algebraic_number_swap(algebraic_roots + to_keep, algebraic_roots + i);
            }
            to_keep ++;
          }
          // Remove the value
          lp_assignment_set_value((lp_assignment_t*) M, x, 0);
          lp_value_destruct(&x_value);
        }
        // Destruct the bad roots
        for (i = to_keep; i < *roots_size; ++ i) {
          lp_algebraic_number_destruct(algebraic_roots + i);
        }
        *roots_size = to_keep;

        // Copy over the roots
        for (i = 0; i < *roots_size; ++i) {
          lp_value_construct(roots + i, LP_VALUE_ALGEBRAIC, algebraic_roots + i);
          lp_algebraic_number_destruct(algebraic_roots + i);
        }
        // Free the temps
        free(algebraic_roots);
        lp_upolynomial_delete(A_alg_u);

      }

      // Remove temps
      coefficient_destruct(&A_alg);

      // Restore the order (remove x from bottom)
      lp_variable_order_make_bot(ctx->var_order, lp_variable_null);
    }
  }

  // Remove temps
  coefficient_destruct(&A_rat);
  integer_destruct(&multiplier);
}

lp_value_t* coefficient_evaluate(const lp_polynomial_context_t* ctx, const coefficient_t* C, const lp_assignment_t* M) {

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_evaluate("); coefficient_print(ctx, C, trace_out); tracef(")\n");
  }

  lp_value_t* result = 0;

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_evaluate(): approximating\n");
  }

  // Compute the initial value of the approximation
  lp_rational_interval_t C_approx;
  lp_rational_interval_construct_zero(&C_approx);
  coefficient_value_approx(ctx, C, M, &C_approx);

  if (trace_is_enabled("coefficient")) {
    tracef("coefficient_evaluate(): approximation = >"); lp_rational_interval_print(&C_approx, trace_out); tracef("\n");
  }

  if (lp_rational_interval_is_point(&C_approx)) {
    // If the approximation is a point we're done
    result = lp_value_new(LP_VALUE_RATIONAL, lp_rational_interval_get_point(&C_approx));
    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_evaluate(): approximation is rational = >"); lp_value_print(result, trace_out); tracef("\n");
    }
  } else {

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_evaluate(): approximation is not rational\n");
    }

    //
    // To evaluate the polynomial C(x1, ..., xn) with values v1, ..., vn
    // we construct the polynomial
    //
    //   A(y, x1, ..., xn) = y - C(x1, ..., xn)
    //
    // and resolve out the defining polynomials p1(x1), ..., pn(xn) to get B(y)
    //
    // B(y) doesn't  vanish since y has coefficient 1.
    //
    // Value of C is one of the roots of B(y). We can isolate the roots and
    // pick the one that intersects with the approximate value of C.

    // Get the variable
    lp_variable_t y = lp_polynomial_context_get_temp_variable(ctx);
    lp_variable_order_make_bot(ctx->var_order, y);

    coefficient_t A;
    lp_integer_t one;
    integer_construct_from_int(lp_Z, &one, 1);
    coefficient_construct_simple(ctx, &A, &one, y, 1);
    coefficient_sub(ctx, &A, &A, C);
    assert(VAR(&A) != y); // y should be the bottom variable

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_evaluate(): resolving algebraic numbers\n");
      tracef("coefficient_evaluate(): start = "); coefficient_print(ctx, &A, trace_out), tracef("\n");
    }

    // Resolve the algebraic numbers
    coefficient_t B;
    coefficient_construct(ctx, &B);
    coefficient_resolve_algebraic(ctx, &A, M, &B);
    assert(coefficient_is_univariate(&B));
    assert(VAR(&B) == y);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_evaluate(): resolved = "); coefficient_print(ctx, &B, trace_out), tracef("\n");
      tracef("coefficient_evaluate(): univariate root isolation\n");
    }

    size_t A_resolved_deg = coefficient_degree(&B);
    lp_value_t* roots = malloc(sizeof(lp_value_t)*A_resolved_deg);
    size_t roots_size = 0;

    // Find the roots of A_resolved
    coefficient_roots_isolate_univariate(ctx, &B, roots, &roots_size);

    if (trace_is_enabled("coefficient")) {
      tracef("coefficient_evaluate(): roots:\n");
      size_t i = 0;
      for (i = 0; i < roots_size; ++ i) {
        tracef("%zu: ", i);
        lp_value_print(roots + i, trace_out);
        tracef("\n");
      }
    }

    // Cache the variables
    lp_variable_list_t C_vars;
    lp_variable_list_construct(&C_vars);
    coefficient_get_variables(C, &C_vars);

    // Cache the intervals
    lp_dyadic_interval_t* rational_interval_cache = algebraic_interval_remember(&C_vars, M);

    // Filter the roots until we get only one
    for (;;) {

      assert(!lp_rational_interval_is_point(&C_approx));

      // Filter
      size_t i, to_keep;
      for (i = 0, to_keep = 0; i < roots_size; ++ i) {
        if (lp_rational_interval_contains_value(&C_approx, roots + i)) {
          lp_value_assign(roots + to_keep, roots + i);
          to_keep ++;
        }
      }
      for (i = to_keep; i < roots_size; ++ i) {
        lp_value_destruct(roots + i);
      }
      roots_size = to_keep;
      assert(roots_size > 0);

      if (trace_is_enabled("coefficient")) {
        tracef("coefficient_evaluate(): filtered roots:\n");
        size_t i = 0;
        for (i = 0; i < roots_size; ++ i) {
          tracef("%zu: ", i);
          lp_value_print(roots + i, trace_out);
          tracef("\n");
        }
      }

      // Did we find it
      if (roots_size == 1) {
        break;
      }

      // Nope, refine
      for (i = 0; i < C_vars.list_size; ++ i) {
        lp_variable_t x_i = C_vars.list[i];
        const lp_value_t* x_i_value = lp_assignment_get_value(M, x_i);
        if (!lp_value_is_rational(x_i_value)) {
          lp_algebraic_number_refine_const(&x_i_value->value.a);
        }
      }

      // Approximate again
      coefficient_value_approx(ctx, C, M, &C_approx);

      // If we reduced to a point, we're done
      if (lp_rational_interval_is_point(&C_approx)) {
        break;
      }
    }

    if (lp_rational_interval_is_point(&C_approx)) {
      // Rational value approximation reduced to
      result = lp_value_new(LP_VALUE_RATIONAL, lp_rational_interval_get_point(&C_approx));
    } else {
      // We have the value in roots[0]
      result = lp_value_new_copy(roots);
    }

    // Remove the roots
    assert(roots_size == 1);
    lp_value_destruct(roots);
    free(roots);

    // Release the variable
    lp_polynomial_context_release_temp_variable(ctx, y);
    lp_variable_order_make_bot(ctx->var_order, lp_variable_null);

    // Release the cache
    algebraic_interval_restore(&C_vars, rational_interval_cache, M);

    // Remove temps
    integer_destruct(&one);
    coefficient_destruct(&A);
    coefficient_destruct(&B);
    lp_variable_list_destruct(&C_vars);

  }

  lp_rational_interval_destruct(&C_approx);

  return result;
}
