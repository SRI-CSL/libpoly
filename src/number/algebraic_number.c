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

#include <interval.h>
#include <upolynomial.h>
#include <algebraic_number.h>
#include <polynomial_context.h>
#include <variable_db.h>

#include "number/integer.h"
#include "interval/arithmetic.h"
#include "polynomial/coefficient.h"
#include "polynomial/output.h"
#include "upolynomial/output.h"

#include "utils/debug_trace.h"

#include <assert.h>

static
void lp_algebraic_number_refine_with_point(const lp_algebraic_number_t* a_const, const lp_dyadic_rational_t* q);

void lp_algebraic_number_construct(lp_algebraic_number_t* a, lp_upolynomial_t* f, const lp_dyadic_interval_t* lr) {
  assert(f);
  // Zero should always constructed separately
  assert(lp_upolynomial_const_term(f));
  assert(lr->a_open && lr->b_open);
  assert(lp_upolynomial_is_primitive(f));

  a->f = f;
  lp_dyadic_interval_construct_copy(&a->I, lr);
  a->sgn_at_a = lp_upolynomial_sgn_at_dyadic_rational(f, &a->I.a);
  a->sgn_at_b = lp_upolynomial_sgn_at_dyadic_rational(f, &a->I.b);
  assert(a->sgn_at_a * a->sgn_at_b < 0);

  // Refine so that the interval is < 1, i.e. log < 0
  while (lp_dyadic_interval_size(&a->I) >= 0) {
    lp_algebraic_number_refine(a);
  }

  // We now have at most one integer in the interval, so just check the options
  if (a->f) {
    lp_dyadic_rational_t value;
    dyadic_rational_construct(&value);
    dyadic_rational_ceiling(&a->I.a, &value);
    lp_algebraic_number_refine_with_point(a, &value);
    dyadic_rational_destruct(&value);
  }
  if (a->f) {
    lp_dyadic_rational_t value;
    dyadic_rational_construct(&value);
    dyadic_rational_floor(&a->I.b, &value);
    lp_algebraic_number_refine_with_point(a, &value);
    dyadic_rational_destruct(&value);
  }
}

void lp_algebraic_number_construct_zero(lp_algebraic_number_t* a) {
  a->f = 0;
  lp_dyadic_interval_construct_from_int(&a->I, 0, 0, 0, 0);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void lp_algebraic_number_construct_one(lp_algebraic_number_t* a) {
  a->f = 0;
  lp_dyadic_interval_construct_from_int(&a->I, 1, 0, 1, 0);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void lp_algebraic_number_construct_copy(lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2) {
  a1->f = a2->f ? lp_upolynomial_construct_copy(a2->f) : 0;
  lp_dyadic_interval_construct_copy(&a1->I, &a2->I);
  a1->sgn_at_a = a2->sgn_at_a;
  a1->sgn_at_b = a2->sgn_at_b;
}

void lp_algebraic_number_construct_from_integer(lp_algebraic_number_t* a, const lp_integer_t* z) {
  lp_dyadic_rational_t z_dy;
  lp_dyadic_rational_construct_from_integer(&z_dy, z);
  lp_algebraic_number_construct_from_dyadic_rational(a, &z_dy);
  lp_dyadic_rational_destruct(&z_dy);
}

void lp_algebraic_number_construct_from_rational(lp_algebraic_number_t* a, const lp_rational_t* rat) {

  // if the rational is integer, we just return it
  if (rational_is_integer(rat)) {
    lp_algebraic_number_construct_from_integer(a, rational_get_num_ref(rat));
    return;
  }

  // x = p/q -> q*x - p = 0
  lp_integer_t c[2];
  lp_integer_t* p_neg = c;
  lp_integer_t* q = c + 1;

  // set the coefficients
  integer_construct(p_neg);
  rational_get_num(rat, p_neg);
  integer_neg(lp_Z, p_neg, p_neg);
  integer_construct(q);
  rational_get_den(rat, q);

  // make the polynomial
  lp_upolynomial_t* f = lp_upolynomial_construct(lp_Z, 1, c);

  // get the floor and ceiling
  lp_integer_t rat_floor, rat_ceil;
  integer_construct(&rat_floor);
  integer_construct(&rat_ceil);
  rational_floor(rat, &rat_floor);
  rational_ceiling(rat, &rat_ceil);

  // since rat is not an integer, we're ok to take these as endpoints
  lp_dyadic_interval_t I;
  lp_dyadic_interval_construct_from_integer(&I, &rat_floor, 1, &rat_ceil, 1);

  // construct the number
  lp_algebraic_number_construct(a, f, &I);

  // remove temps
  lp_dyadic_interval_destruct(&I);
  integer_destruct(&rat_floor);
  integer_destruct(&rat_ceil);
  integer_destruct(q);
  integer_destruct(p_neg);
}


void lp_algebraic_number_construct_from_dyadic_rational(lp_algebraic_number_t* a, const lp_dyadic_rational_t* q) {
  a->f = 0;
  lp_dyadic_interval_construct_point(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

void lp_algebraic_number_destruct(lp_algebraic_number_t* a) {
  if (a->f) {
    lp_upolynomial_delete(a->f);
  }
  lp_dyadic_interval_destruct(&a->I);
}

void lp_algebraic_number_swap(lp_algebraic_number_t* a, lp_algebraic_number_t* b) {
  lp_algebraic_number_t tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline
void lp_algebraic_number_reduce_polynomial(const lp_algebraic_number_t* a, const lp_upolynomial_t* f, int sgn_at_a, int sgn_at_b) {
  __var_unused(sgn_at_a);
  __var_unused(sgn_at_b);
  assert(a->f);
  assert(a->sgn_at_a * a->sgn_at_b < 0);
  assert(sgn_at_a * sgn_at_b < 0);
  assert(lp_upolynomial_is_primitive(f));
  lp_algebraic_number_t* a_nonconst = (lp_algebraic_number_t*) a;
  lp_upolynomial_delete(a_nonconst->f);
  a_nonconst->f = lp_upolynomial_construct_copy(f);
  a_nonconst->sgn_at_a = sgn_at_a;
  a_nonconst->sgn_at_b = sgn_at_b;
}

static inline
void lp_algebraic_number_collapse_to_point(const lp_algebraic_number_t* a_const, const lp_dyadic_rational_t* q) {
  assert(a_const->f);
  assert(lp_upolynomial_sgn_at_dyadic_rational(a_const->f, q) == 0);
  // We'll modify the number so unconst it
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  lp_upolynomial_delete(a->f);
  a->f = 0;
  lp_dyadic_interval_collapse_to(&a->I, q);
  a->sgn_at_a = 0;
  a->sgn_at_b = 0;
}

static inline
void lp_algebraic_number_refine_with_point(const lp_algebraic_number_t* a_const, const lp_dyadic_rational_t* q) {
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  if (a->f && lp_dyadic_interval_contains_dyadic_rational(&a->I, q)) {
    // Compute the sign at the left end point
    int a_sgn_at_q = lp_upolynomial_sgn_at_dyadic_rational(a->f, q);
    if (a_sgn_at_q == 0) {
      lp_algebraic_number_collapse_to_point(a, q);
      return;
    } else if (a_sgn_at_q * a->sgn_at_a > 0) {
      // We can keep I.a for a's a
      lp_dyadic_interval_set_a(&a->I, q, 1);
    } else  {
      // We can keep I.a for a's b
      lp_dyadic_interval_set_b(&a->I, q, 1);
    }
  }
}

/**
 * Refine the interval using one split. Returns -1 if left, +1 if right, 0 if
 * reduced to point.
 */
static inline
int lp_algebraic_number_refine_const_internal(const lp_algebraic_number_t* a_const) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine(");
    lp_algebraic_number_print(a_const, trace_out);
    tracef(")\n");
  }

  assert(a_const->f);

  int result;

  // We'll modify the number so unconst it
  lp_algebraic_number_t* a = (lp_algebraic_number_t*) a_const;
  // Compute the mid point
  lp_dyadic_interval_t I_left;
  lp_dyadic_interval_t I_right;
  lp_dyadic_interval_construct_from_split(&I_left, &I_right, &a_const->I, 1, 1);
  // Compute the sign at the mid-point
  const lp_dyadic_rational_t* m = &I_left.b;
  int sgn_at_m = lp_upolynomial_sgn_at_dyadic_rational(a_const->f, m);
  if (sgn_at_m == 0) {
    // m is actually the number a1
    lp_algebraic_number_collapse_to_point(a_const, m);
    result = 0;
  } else if (sgn_at_m * a_const->sgn_at_a > 0) {
    // a m b  or  a m b
    // + + -      - - +
    lp_dyadic_interval_swap(&I_right, &a->I);
    result = 1;
  } else {
    // a m b  or  a m b
    // + - -      - + +
    lp_dyadic_interval_swap(&I_left, &a->I);
    result = -1;
  }
  // Remove temp
  lp_dyadic_interval_destruct(&I_left);
  lp_dyadic_interval_destruct(&I_right);

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_refine() => ");
    lp_algebraic_number_print(a_const, trace_out);
    tracef(", d = %d\n", result);
  }

  return result;
}

void lp_algebraic_number_refine(lp_algebraic_number_t* a) {
  if (a->f) {
    lp_algebraic_number_refine_const_internal(a);
  }
}

void lp_algebraic_number_refine_const(const lp_algebraic_number_t* a) {
  if (a->f) {
    lp_algebraic_number_refine_const_internal(a);
  }
}

void lp_algebraic_number_restore_interval(lp_algebraic_number_t* a, const lp_dyadic_interval_t* I) {
  lp_dyadic_interval_assign(&a->I, I);
}

void lp_algebraic_number_restore_interval_const(const lp_algebraic_number_t* a_const, const lp_dyadic_interval_t* I) {
  lp_algebraic_number_restore_interval((lp_algebraic_number_t* ) a_const, I);
}


int lp_algebraic_number_sgn(const lp_algebraic_number_t* a) {
  lp_integer_t zero;
  lp_integer_construct_from_int(lp_Z, &zero, 0);
  int sgn = lp_algebraic_number_cmp_integer(a, &zero);
  lp_integer_destruct(&zero);
  return sgn;
}

/**
 * The "proper" algebraic number is always a1.
 */
int lp_algebraic_number_cmp(const lp_algebraic_number_t* a1, const lp_algebraic_number_t* a2) {

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // We only have a problem if the intervals intersect
  if (!lp_dyadic_interval_disjoint(&a1->I, &a2->I)) {
    // First intersect the intervals. Since both intervals are open or points
    // the intersection is either open or a point
    lp_dyadic_interval_t I;
    lp_dyadic_interval_construct_intersection(&I, &a1->I, &a2->I);

    if (trace_is_enabled("algebraic_number")) {
      tracef("I = "); lp_dyadic_interval_print(&I, trace_out); tracef("\n");
    }

    // Refine the interval using the intersection points
    lp_algebraic_number_refine_with_point(a1, &I.a);
    lp_algebraic_number_refine_with_point(a2, &I.a);
    if (!I.is_point) {
      lp_algebraic_number_refine_with_point(a1, &I.b);
      lp_algebraic_number_refine_with_point(a2, &I.b);
    }

    // Remove the temp intersection
    lp_dyadic_interval_destruct(&I);
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(")\n");
  }

  // At this point the intervals can intersect only if they are equal
  // In this case
  int equal = 0;
  if (a1->f && a2->f && lp_dyadic_interval_equals(&a1->I, &a2->I)) {
    // If we are of the same size, and equal intervals, check for equality
    // We check by first intersecting the two
    // Both proper algebraic, intervals are equal
    // Let f = f'*gcd, g = g'*gcd
    // Gcd has a zero in the interval iff the numbers are equal and
    // we can reduce polynomials to the gcd
    lp_upolynomial_t* gcd = lp_upolynomial_gcd(a1->f, a2->f);
    int sgn_at_a = lp_upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.a);
    int sgn_at_b = lp_upolynomial_sgn_at_dyadic_rational(gcd, &a1->I.b);
    if (sgn_at_a * sgn_at_b < 0) {
      lp_algebraic_number_reduce_polynomial(a1, gcd, sgn_at_a, sgn_at_b);
      lp_algebraic_number_reduce_polynomial(a2, gcd, sgn_at_a, sgn_at_b);
      equal = 1;
    } else {
      // We're not equal, so bisect away
      int d1 = 1, d2 = 1;
      while (d1 == d2 && d1 && d2) {
        // They become different when bisection goes different ways
        d1 = lp_algebraic_number_refine_const_internal(a1);
        d2 = lp_algebraic_number_refine_const_internal(a2);
      }
    }
    lp_upolynomial_delete(gcd);
  }

  int result;

  if (equal) {
    // We checked for equality and they are equal
    result = 0;
  } else {
    // Not equal, but disjunct, so compare the two intervals
    int cmp = dyadic_rational_cmp(&a1->I.a, &a2->I.a);
    if (cmp == 0) {
      // Since they are disjunct one of them is open
      if (a1->I.a_open && !a2->I.a_open) {
        // a1  (     )
        // a2  []
        result = 1;
      } else if (!a1->I.a_open && a2->I.a_open) {
        // a1  []
        // a2  (    )
        result = -1;
      } else {
        result = cmp;
      }
    } else {
      result = cmp;
    }
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("algebraic_number_cmp(");
    lp_algebraic_number_print(a1, trace_out);
    tracef(", ");
    lp_algebraic_number_print(a2, trace_out);
    tracef(") => %d\n", result);
  }

  return result;
}

int lp_algebraic_number_cmp_void(const void* a1, const void* a2) {
  return lp_algebraic_number_cmp(a1, a2);
}

int lp_algebraic_number_cmp_integer(const lp_algebraic_number_t* a1, const lp_integer_t* a2) {
  if (a1->f) {
    assert(!a1->I.is_point);
    // Easy check, compare to the dyadic interval
    int cmp = lp_dyadic_interval_cmp_integer(&a1->I, a2);
    if (cmp != 0) {
      return cmp;
    }
    // Point in interval, let's evaluate
    int poly_sgn = lp_upolynomial_sgn_at_integer(a1->f, a2);
    if (poly_sgn == 0) {
      return 0;
    }
    // Not a zero, so bisect while not outside
    while (cmp == 0) {
      lp_algebraic_number_refine_const_internal(a1);
      cmp = lp_dyadic_interval_cmp_integer(&a1->I, a2);
    }
    // Return the last compare
    return cmp;
  } else {
    assert(a1->I.is_point);
    return dyadic_rational_cmp_integer(&a1->I.a, a2);
  }
}

int lp_algebraic_number_cmp_dyadic_rational(const lp_algebraic_number_t* a1, const lp_dyadic_rational_t* a2) {
  if (a1->f) {
    assert(!a1->I.is_point);
    // Easy check, compare to the dyadic interval
    int cmp = lp_dyadic_interval_cmp_dyadic_rational(&a1->I, a2);
    if (cmp != 0) {
      return cmp;
    }
    // Point in interval, let's evaluate
    int poly_sgn = lp_upolynomial_sgn_at_dyadic_rational(a1->f, a2);
    if (poly_sgn == 0) {
      return 0;
    }
    // Not a zero, so bisect while not outside
    while (cmp == 0) {
      lp_algebraic_number_refine_const_internal(a1);
      cmp = lp_dyadic_interval_cmp_dyadic_rational(&a1->I, a2);
    }
    // Return the last compare
    return cmp;
  } else {
    assert(a1->I.is_point);
    return dyadic_rational_cmp(&a1->I.a, a2);
  }
}

int lp_algebraic_number_cmp_rational(const lp_algebraic_number_t* a1, const lp_rational_t* a2) {
  if (a1->f) {
    assert(!a1->I.is_point);
    // Easy check, compare to the dyadic interval
    int cmp = lp_dyadic_interval_cmp_rational(&a1->I, a2);
    if (cmp != 0) {
      return cmp;
    }
    // Point in interval, let's evaluate
    int poly_sgn = lp_upolynomial_sgn_at_rational(a1->f, a2);
    if (poly_sgn == 0) {
      return 0;
    }
    // Not a zero, so bisect while not outside
    while (cmp == 0) {
      lp_algebraic_number_refine_const_internal(a1);
      cmp = lp_dyadic_interval_cmp_rational(&a1->I, a2);
    }
    // Return the last compare
    return cmp;
  } else {
    return -rational_cmp_dyadic_rational(a2, &a1->I.a);
  }
}


int lp_algebraic_number_print(const lp_algebraic_number_t* a, FILE* out) {
  if (a->f == 0) {
    return dyadic_rational_print(&a->I.a, out);
  } else {
    int ret = 0;
    ret += fprintf(out, "<");
    ret += lp_upolynomial_print(a->f, out);
    ret += fprintf(out, ", ");
    ret += lp_dyadic_interval_print(&a->I, out);
    ret += fprintf(out, ">");
    return ret;
  }
}

char* lp_algebraic_number_to_string(const lp_algebraic_number_t* a) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_algebraic_number_print(a, f);
  fclose(f);
  return str;
}

double lp_algebraic_number_to_double(const lp_algebraic_number_t* a_const) {

  // If a point, just return it's double
  if (a_const->f == 0) {
    return dyadic_rational_to_double(&a_const->I.a);
  }

  // We do the necessary refinement on a copy
  lp_algebraic_number_t a;
  lp_algebraic_number_construct_copy(&a, a_const);

  // Refine the number until we get the desired precision
  lp_dyadic_rational_t interval_size;
  dyadic_rational_construct(&interval_size);
  dyadic_rational_sub(&interval_size, &a.I.b, &a.I.a);
  if (interval_size.n < 100) {
    int iterations = 100 - interval_size.n;
    while (a.f && iterations > 0) {
      lp_algebraic_number_refine_const_internal(&a);
      iterations --;
    }
  }

  double result = dyadic_rational_to_double(&a.I.a);

  dyadic_rational_destruct(&interval_size);
  lp_algebraic_number_destruct(&a);

  return result;
}

void lp_algebraic_number_to_rational(const lp_algebraic_number_t* a_const, lp_rational_t* q) {

  lp_rational_t tmp;

  // If a point, just return it's double
  if (a_const->f == 0) {
    rational_construct_from_dyadic(&tmp, &a_const->I.a);
    rational_swap(q, &tmp);
    rational_destruct(&tmp);
    return;
  }

  // We do the necessary refinement on a copy
  lp_algebraic_number_t a;
  lp_algebraic_number_construct_copy(&a, a_const);

  // Refine the number until we get the desired precision
  lp_dyadic_rational_t interval_size;
  dyadic_rational_construct(&interval_size);
  dyadic_rational_sub(&interval_size, &a.I.b, &a.I.a);
  if (interval_size.n < 100) {
    int iterations = 100 - interval_size.n;
    while (a.f && iterations > 0) {
      lp_algebraic_number_refine_const_internal(&a);
      iterations --;
    }
  }

  rational_construct_from_dyadic(&tmp, &a.I.a);
  rational_swap(q, &tmp);
  rational_destruct(&tmp);

  dyadic_rational_destruct(&interval_size);
  lp_algebraic_number_destruct(&a);
}


/** Polynomial context for resultant computation */
static const lp_polynomial_context_t* algebraic_ctx = 0;

/** Result variable */
static lp_variable_t var_r;

/** Variable for left operand */
static lp_variable_t var_x;

/** Variable for right operand */
static lp_variable_t var_y;

static
const lp_polynomial_context_t* lp_algebraic_pctx(void) {
  // Create the context if needed
  if (algebraic_ctx == 0) {
    // Create the variables
    lp_variable_db_t* var_db = lp_variable_db_new();
    var_x = lp_variable_db_new_variable(var_db, "_x");
    var_y = lp_variable_db_new_variable(var_db, "_y");
    var_r = lp_variable_db_new_variable(var_db, "_r");
    // Order as Z[r, y, x]
    lp_variable_order_t* var_order = lp_variable_order_new();
    lp_variable_order_push(var_order, var_r);
    lp_variable_order_push(var_order, var_y);
    lp_variable_order_push(var_order, var_x);
    // Create the context
    algebraic_ctx = lp_polynomial_context_new(0, var_db, var_order);
    // Detach local references
    lp_variable_db_detach(var_db);
    lp_variable_order_detach(var_order);
  }
  return algebraic_ctx;
}

void filter_roots(lp_algebraic_number_t* roots, size_t* roots_size, const lp_dyadic_interval_t* I) {
  size_t i, to_keep;
  for (i = 0, to_keep = 0; i < *roots_size; ++ i) {
    lp_algebraic_number_t* root = roots + i;
    if (lp_dyadic_interval_disjoint(&root->I, I)) {
      // Remove this root
      lp_algebraic_number_destruct(root);
    } else {
      // Keep it
      if (i > to_keep) {
        *(roots + to_keep) = *root;
      }
      to_keep ++;
    }
  }
  *roots_size = to_keep;
}

/** Function type called on coefficient traversal (such as r - (x + y)) */
typedef void (*construct_op_polynomial_f) (coefficient_t* op, void* data);

/** Function type called on interval operations (such as I = I1 + I2) */
typedef void (*interval_op_f) (lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data);

/** Compute the operation */
static
void lp_algebraic_number_op(
    lp_algebraic_number_t* op, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b,
    construct_op_polynomial_f construct_op,
    interval_op_f interval_op,
    void* data)
{
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  if (trace_is_enabled("algebraic_number")) {
    tracef("a = "); lp_algebraic_number_print(a, trace_out); tracef("\n");
    if (b) {
      tracef("b = "); lp_algebraic_number_print(b, trace_out); tracef("\n");
    }
    tracef("op = "); lp_algebraic_number_print(op, trace_out); tracef("\n");
  }

  coefficient_t f_a;
  if (a->f) {
    coefficient_construct_from_univariate(ctx, &f_a, a->f, var_x);
  } else {
    assert(a->I.is_point);
    // x = p/q -> q*x - p = 0
    lp_integer_t p_neg, q;
    integer_construct(&p_neg);
    integer_construct(&q);
    integer_neg(lp_Z, &p_neg, &a->I.a.a);
    dyadic_rational_get_den(&a->I.a, &q);
    coefficient_construct_linear(ctx, &f_a, &q, &p_neg, var_x);
    lp_integer_destruct(&p_neg);
    lp_integer_destruct(&q);
  }
  if (trace_is_enabled("algebraic_number")) {
    tracef("f_a = "); coefficient_print(ctx, &f_a, trace_out); tracef("\n");
  }

  coefficient_t f_b;
  if (b) {
    if (b->f) {
      coefficient_construct_from_univariate(ctx, &f_b, b->f, var_y);
    } else {
      // x = p/q -> q*x - p = 0
      lp_integer_t p_neg, q;
      integer_construct(&p_neg);
      integer_construct(&q);
      integer_neg(lp_Z, &p_neg, &b->I.a.a);
      dyadic_rational_get_den(&b->I.a, &q);
      coefficient_construct_linear(ctx, &f_b, &q, &p_neg, var_y);
      lp_integer_destruct(&p_neg);
      lp_integer_destruct(&q);
    }
    if (trace_is_enabled("algebraic_number")) {
      tracef("f_b = "); coefficient_print(ctx, &f_b, trace_out); tracef("\n");
    }
  }

  // Construct the op polynomial
  coefficient_t f_r;
  construct_op(&f_r, data);

  if (trace_is_enabled("algebraic_number")) {
    tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
  }

  // Compute the resultant
  coefficient_resultant(ctx, &f_r, &f_r, &f_a);
  if (trace_is_enabled("algebraic_number")) {
    tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
  }
  if (b) {
    coefficient_resultant(ctx, &f_r, &f_r, &f_b);
    if (trace_is_enabled("algebraic_number")) {
      tracef("f_r = "); coefficient_print(ctx, &f_r, trace_out); tracef("\n");
    }
  }

  // Resultant polynomial captures the result
  lp_upolynomial_t* f = coefficient_to_univariate(ctx, &f_r);

  if (trace_is_enabled("algebraic_number")) {
    tracef("f = "); lp_upolynomial_print(f, trace_out);
  }

  // Get the roots of f
  size_t f_roots_size = 0;
  lp_algebraic_number_t* f_roots = malloc(sizeof(lp_algebraic_number_t)*lp_upolynomial_degree(f));
  lp_upolynomial_roots_isolate(f, f_roots, &f_roots_size);
  lp_upolynomial_delete(f);

  // Interval for the result
  lp_dyadic_interval_t I;
  lp_dyadic_interval_construct_zero(&I);

  // Approximate the result
  while (f_roots_size > 1) {

    // Add the two intervals
    if (b) {
      interval_op(&I, &a->I, &b->I, data);
    } else {
      interval_op(&I, &a->I, 0, data);
    }

    if (trace_is_enabled("algebraic_number")) {
      tracef("a = "); lp_algebraic_number_print(a, trace_out); tracef("\n");
      if (b) {
        tracef("b = "); lp_algebraic_number_print(b, trace_out); tracef("\n");
      }
      tracef("I = "); lp_dyadic_interval_print(&I, trace_out); tracef("\n");
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        tracef("f[%zu] = ", i); lp_algebraic_number_print(f_roots + i, trace_out); tracef("\n");
      }
    }

    // Filter the roots over I
    filter_roots(f_roots, &f_roots_size, &I);

    // If more then one root, we need to refine a and b
    if (f_roots_size > 1) {
      if (a->f) {
        lp_algebraic_number_refine_const_internal(a);
      }
      if (b && b->f) {
        lp_algebraic_number_refine_const_internal(b);
      }
      size_t i;
      for (i = 0; i < f_roots_size; ++ i) {
        if ((f_roots + i)->f) {
          // not a point, refine it
          lp_algebraic_number_refine_const_internal(f_roots + i);
        }
      }
    }
  }

  assert(f_roots_size == 1);

  // Our root is the only one left
  lp_algebraic_number_destruct(op);
  *op = *f_roots;

  if (trace_is_enabled("algebraic_number")) {
    tracef("op = "); lp_algebraic_number_print(op, trace_out); tracef("\n");
  }

  // Remove temps
  coefficient_destruct(&f_a);
  if (b) {
    coefficient_destruct(&f_b);
  }
  coefficient_destruct(&f_r);
  lp_dyadic_interval_destruct(&I);
  free(f_roots);
}

static
void lp_algebraic_number_add_construct_op(coefficient_t* f_r, void* data) {
  __var_unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x + y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  coefficient_sub(ctx, f_r, f_r, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

static
void lp_algebraic_number_add_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __var_unused(data);
  dyadic_interval_add(I, I1, I2);
}

lp_upolynomial_t* lp_upolynomial_shift(const lp_upolynomial_t* poly, const lp_integer_t* shift_num, unsigned long shift_den_exp) {
  lp_integer_t shift_den;
  lp_integer_construct_from_int(lp_Z, &shift_den, -1);
  lp_integer_mul_pow2(lp_Z, &shift_den, &shift_den, shift_den_exp);
  lp_integer_t multcoeffs[2] = { *shift_num, shift_den };
  lp_upolynomial_t* basemult = lp_upolynomial_construct(lp_Z, 1, multcoeffs);
  lp_upolynomial_t* curmult = lp_upolynomial_construct_copy(basemult);

  lp_integer_t coeffs[lp_upolynomial_degree(poly) + 1];
  for (size_t i = 0; i <= lp_upolynomial_degree(poly); ++i) {
    lp_integer_construct(&coeffs[i]);
  }
  lp_upolynomial_unpack(poly, coeffs);

  lp_upolynomial_t* out = lp_upolynomial_construct(lp_Z, 0, coeffs);
  
  for (size_t i = 1; i <= lp_upolynomial_degree(poly); ++i) {
    lp_upolynomial_t* cur = lp_upolynomial_mul_c(curmult, &coeffs[i]);
    
    lp_upolynomial_t* tout = lp_upolynomial_add(out, cur);
    lp_upolynomial_t* tcurmult = lp_upolynomial_mul(curmult, basemult);
    lp_upolynomial_delete(cur);
    lp_upolynomial_delete(curmult);
    lp_upolynomial_delete(out);
    out = tout;
    curmult = tcurmult;
  }
  for (size_t i = 0; i <= lp_upolynomial_degree(poly); ++i) {
    lp_integer_destruct(&coeffs[i]);
  }
  lp_integer_destruct(&shift_den);

  lp_upolynomial_delete(basemult);
  lp_upolynomial_delete(curmult);
  return out;
}

void lp_algebraic_number_add(lp_algebraic_number_t* sum, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  lp_algebraic_number_op(sum, a, b, lp_algebraic_number_add_construct_op, lp_algebraic_number_add_interval_op, 0);
}

void lp_algebraic_number_sub_construct_op(coefficient_t* f_r, void* data) {
  __var_unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x - y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  coefficient_add(ctx, f_r, f_r, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

void lp_algebraic_number_sub_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __var_unused(data);
  dyadic_interval_sub(I, I1, I2);
}

void lp_algebraic_number_sub(lp_algebraic_number_t* sub, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  lp_algebraic_number_op(sub, a, b, lp_algebraic_number_sub_construct_op, lp_algebraic_number_sub_interval_op, 0);
}

void lp_algebraic_number_neg(lp_algebraic_number_t* neg, const lp_algebraic_number_t* a) {
  lp_upolynomial_t* f_neg_x = 0;

  if (a->f) {
    f_neg_x = lp_upolynomial_subst_x_neg(a->f);
    if (integer_sgn(lp_Z, lp_upolynomial_lead_coeff(f_neg_x)) < 0) {
      lp_upolynomial_neg_in_place(f_neg_x);
    }
  }

  lp_dyadic_interval_t I_neg; // To destroy
  lp_dyadic_interval_construct_copy(&I_neg, &a->I);
  dyadic_interval_neg(&I_neg, &I_neg);

  lp_algebraic_number_t result; // To destroy
  lp_algebraic_number_construct(&result, f_neg_x, &I_neg);
  lp_algebraic_number_swap(&result, neg);

  lp_algebraic_number_destruct(&result);
  lp_dyadic_interval_destruct(&I_neg);
}

void lp_algebraic_number_mul_construct_op(coefficient_t* f_r, void* data) {
  __var_unused(data);
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  // Construct the polynomial z - (x*y)
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x, f_y;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, 1);
  coefficient_construct_simple(ctx, &f_y, &one, var_y, 1);
  coefficient_sub_mul(ctx, f_r, &f_x, &f_y);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
  coefficient_destruct(&f_y);
}

void lp_algebraic_number_mul_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  __var_unused(data);
  dyadic_interval_mul(I, I1, I2);
}

void lp_algebraic_number_mul(lp_algebraic_number_t* mul, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  // Special case, when one is zero
  if (lp_algebraic_number_sgn(a) == 0 || lp_algebraic_number_sgn(b) == 0) {
    lp_algebraic_number_destruct(mul);
    lp_algebraic_number_construct_zero(mul);
  } else {
    lp_algebraic_number_op(mul, a, b, lp_algebraic_number_mul_construct_op, lp_algebraic_number_mul_interval_op, 0);
  }
}

/**
 * Multiplicative inverse (a != 0)
 *
 * f(x) = a_n x^n + ... + a_0
 *
 * with f(a) = 0, a only root in (l, u), 0 not in (l, u)
 *
 * NOTE: l, u could be 0
 *
 * g(x) = x^n*f(1/x) = x^n*(a_n/x^n + ... + a_0)
 *      = a_n + ... + a_0 x^n
 *
 * with g(1/a) = 0, only root in 1/u < 1/a < 1/l
 *
 * but, 1/u, 1/l are not dyadic rationals...
 *
 * g(x) is square-free otherwise g(x) = h1(x)^2*h2(x) and so
 * f(x) = x^n*g(1/x) would not be square free
 */
void lp_algebraic_number_inv(lp_algebraic_number_t* inv, const lp_algebraic_number_t* a) {

  assert(lp_algebraic_number_sgn(a) != 0);

  if (trace_is_enabled("algebraic_number")) {
    tracef("a = "); lp_algebraic_number_print(a, trace_out); tracef("\n");
  }

  if (a->f == 0) {
    // inverse of dyadic rational is not necessarily dyadic, we just
    // invert the rational and construct from there
    lp_rational_t a_inv_q;
    lp_algebraic_number_t a_inv;
    lp_rational_construct_from_dyadic(&a_inv_q, &a->I.a);
    rational_inv(&a_inv_q, &a_inv_q);
    lp_algebraic_number_construct_from_rational(&a_inv, &a_inv_q);

    // store result
    lp_algebraic_number_swap(inv, &a_inv);

    // remove temps
    lp_algebraic_number_destruct(&a_inv);
    lp_rational_destruct(&a_inv_q);
  } else {

    // refine the number until 0 is not an endpoint of the interval
    while (lp_dyadic_rational_sgn(&a->I.a) == 0 || lp_dyadic_rational_sgn(&a->I.b) == 0) {
      lp_algebraic_number_refine_const(a);
      if (a->f == 0) {
        lp_algebraic_number_inv(inv, a);
        return;
      }
    }

    // construct the inverse polynomial
    lp_upolynomial_t* f_a_inv = lp_upolynomial_construct_copy(a->f);
    lp_upolynomial_reverse_in_place(f_a_inv);
    // Make the polynomial primitive
    if (integer_sgn(lp_Z, lp_upolynomial_lead_coeff(f_a_inv)) < 0) {
      lp_upolynomial_neg_in_place(f_a_inv);
    }

    // we need to construct the interval within (1/u, 1/l) and refine it
    // until we find dyadic end-points sufficient to represent the number
    lp_dyadic_interval_t I;
    lp_dyadic_rational_t I_lb_dy, I_ub_dy;

    // we find sub-intervals in (1/u, 1/l) that are of constant sign
    int sgn;
    lp_rational_t m, I_lb, I_ub;
    rational_construct(&m);
    rational_construct_from_dyadic(&I_lb, &a->I.b);
    rational_construct_from_dyadic(&I_ub, &a->I.a);
    rational_inv(&I_lb, &I_lb);
    rational_inv(&I_ub, &I_ub);

    // find subinterval (I_lb, m) of the same sign as I_lb
    sgn = lp_upolynomial_sgn_at_rational(f_a_inv, &I_lb);
    assert(sgn != 0);
    rational_assign(&m, &I_ub);
    do {
      rational_add(&m, &I_lb, &m);
      rational_div_2exp(&m, &m, 1);
    } while (lp_upolynomial_sgn_at_rational(f_a_inv, &m) * sgn <= 0);
    // pick a dyadic in (I_lb, m)
    dyadic_rational_construct(&I_lb_dy);
    dyadic_rational_get_value_between(&I_lb_dy, &I_lb, &m);

    // find subinterval (m, I_ub) of the same sign as I_ub
    sgn = lp_upolynomial_sgn_at_rational(f_a_inv, &I_ub);
    assert(sgn != 0);
    rational_assign(&m, &I_lb);
    do {
      rational_add(&m, &m, &I_ub);
      rational_div_2exp(&m, &m, 1);
    } while (lp_upolynomial_sgn_at_rational(f_a_inv, &m) * sgn <= 0);
    // oick a dyadic in (m, I_ub)
    dyadic_rational_construct(&I_ub_dy);
    dyadic_rational_get_value_between(&I_ub_dy, &m, &I_ub);

    // construct the root interval
    lp_dyadic_interval_construct(&I, &I_lb_dy, 1, &I_ub_dy, 1);

    // construct the result
    lp_algebraic_number_t a_inv;
    lp_algebraic_number_construct(&a_inv, f_a_inv, &I);

    // store result
    lp_algebraic_number_swap(&a_inv, inv);

    // remove temps
    lp_dyadic_interval_destruct(&I);
    dyadic_rational_destruct(&I_ub_dy);
    dyadic_rational_destruct(&I_lb_dy);
    rational_destruct(&m);
    rational_destruct(&I_lb);
    rational_destruct(&I_ub);
    lp_algebraic_number_destruct(&a_inv);
  }

  if (trace_is_enabled("algebraic_number")) {
    tracef("inv = "); lp_algebraic_number_print(inv, trace_out); tracef("\n");
  }
}

void lp_algebraic_number_div(lp_algebraic_number_t* div, const lp_algebraic_number_t* a, const lp_algebraic_number_t* b) {
  lp_algebraic_number_t inv;
  lp_algebraic_number_construct_zero(&inv);
  lp_algebraic_number_inv(&inv, b);
  lp_algebraic_number_mul(div, a, &inv);
  lp_algebraic_number_destruct(&inv);
}


void lp_algebraic_number_pow_construct_op(coefficient_t* f_r, void* data) {
  const lp_polynomial_context_t* ctx = lp_algebraic_pctx();

  unsigned n = *((unsigned*) data);

  // Construct the polynomial z - x**n
  lp_integer_t one;
  integer_construct_from_int(lp_Z, &one, 1);
  coefficient_t f_x;
  coefficient_construct_simple(ctx, f_r, &one, var_r, 1);
  coefficient_construct_simple(ctx, &f_x, &one, var_x, n);
  coefficient_sub(ctx, f_r, f_r, &f_x);
  integer_destruct(&one);
  coefficient_destruct(&f_x);
}

void lp_algebraic_number_pow_interval_op(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2, void* data) {
  assert(I2 == 0);
  __var_unused(I2);
  unsigned n = *((unsigned*) data);
  dyadic_interval_pow(I, I1, n);
}

void lp_algebraic_number_pow(lp_algebraic_number_t* pow, const lp_algebraic_number_t* a, unsigned n) {
  if (n == 0) {
    // special case x^0 == 1
    lp_integer_t one;
    lp_algebraic_number_t result;
    lp_integer_construct_from_int(lp_Z, &one, 1);
    lp_algebraic_number_construct_from_integer(&result, &one);
    lp_algebraic_number_swap(pow, &result);
    lp_algebraic_number_destruct(&result);
    lp_integer_destruct(&one);
  } else {
    lp_algebraic_number_op(pow, a, 0, lp_algebraic_number_pow_construct_op, lp_algebraic_number_pow_interval_op, &n);
  }
}

void lp_algebraic_number_get_dyadic_midpoint(const lp_algebraic_number_t* a, lp_dyadic_rational_t* q) {
  if (a->I.is_point) {
    lp_dyadic_rational_assign(q, &a->I.a);
  } else {
    lp_dyadic_rational_add(q, &a->I.a, &a->I.b);
    lp_dyadic_rational_div_2exp(q, q, 1);
  }
}

void lp_algebraic_number_get_rational_midpoint(const lp_algebraic_number_t* a, lp_rational_t* q) {
  // Get the dyadic midpoint
  lp_dyadic_rational_t tmp_dy;
  lp_dyadic_rational_construct(&tmp_dy);
  lp_algebraic_number_get_dyadic_midpoint(a, &tmp_dy);
  lp_rational_t tmp_q;
  // Now, convert to rational
  lp_rational_construct_from_dyadic(&tmp_q, &tmp_dy);
  lp_rational_swap(&tmp_q, q);
  lp_rational_destruct(&tmp_q);
  lp_dyadic_rational_destruct(&tmp_dy);
}

int lp_algebraic_number_is_rational(const lp_algebraic_number_t* a) {
  if (lp_dyadic_interval_is_point(&a->I)) {
    // If a point, we're (dyadic) rational
    return 1;
  } else if (lp_upolynomial_degree(a->f) == 1) {
    // If degree 1, we're directly rational
    return 1;
  } else {
    return 0;
  }
}

int lp_algebraic_number_is_integer(const lp_algebraic_number_t* a) {
  if (lp_dyadic_interval_is_point(&a->I)) {
    // If a point, we're (dyadic) rational
    return lp_dyadic_rational_is_integer(&a->I.a);
  } else {
    return 0;
  }
}

void lp_algebraic_number_ceiling(const lp_algebraic_number_t* a, lp_integer_t* a_ceiling) {
  if (lp_dyadic_interval_is_point(&a->I)) {
    dyadic_rational_ceiling_int(&a->I.a, a_ceiling);
  } else {
    dyadic_rational_ceiling_int(&a->I.b, a_ceiling);
  }
}

void lp_algebraic_number_floor(const lp_algebraic_number_t* a, lp_integer_t* a_floor) {
  dyadic_rational_floor_int(&a->I.a, a_floor);
}

size_t lp_algebraic_number_hash_approx(const lp_algebraic_number_t* a, unsigned precision) {

  if (lp_algebraic_number_is_integer(a)) {
    return integer_hash(&a->I.a.a);
  }

  unsigned i;
  size_t hash;
  lp_integer_t a_floor, a_ceil;
  lp_dyadic_rational_t lb, m, ub;

  integer_construct(&a_floor);
  integer_construct(&a_ceil);

  lp_algebraic_number_floor(a, &a_floor);
  lp_algebraic_number_ceiling(a, &a_ceil);

  lp_dyadic_rational_construct_from_integer(&lb, &a_floor);
  lp_dyadic_rational_construct_from_integer(&ub, &a_ceil);
  lp_dyadic_rational_construct_from_integer(&m, &a_floor);

  // refine to precision or until a dyadic rational is found
  for (i = 0; i < precision; ++ i) {
    // m = (lb + ub) / 2
    lp_dyadic_rational_add(&m, &lb, &ub);
    lp_dyadic_rational_div_2exp(&m, &m, 1);
    int cmp = lp_algebraic_number_cmp_dyadic_rational(a, &m);
    if (cmp == 0) {
      break;
    } else if (cmp < 0) {
      // lb < a < m < ub, set m = ub
      lp_dyadic_rational_swap(&m, &ub);
    } else {
      // lb < m < a < ub, set m = lb
      lp_dyadic_rational_swap(&m, &lb);
    }
  }

  // just hash the approximation
  hash = lp_dyadic_rational_hash(&m);

  // remove temps
  lp_dyadic_rational_destruct(&m);
  lp_dyadic_rational_destruct(&ub);
  lp_dyadic_rational_destruct(&lb);
  integer_destruct(&a_ceil);
  integer_destruct(&a_floor);

  return hash;
}
