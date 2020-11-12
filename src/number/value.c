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

#include <rational_interval.h>
#include <algebraic_number.h>
#include <upolynomial.h>

#include "number/value.h"
#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

#include "utils/debug_trace.h"

#include <limits.h>
#include <math.h>

void lp_value_construct(lp_value_t* v, lp_value_type_t type, const void* data) {
  v->type = type;
  switch(type) {
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
  case LP_VALUE_NONE:
    break;
  case LP_VALUE_INTEGER:
    integer_construct_copy(lp_Z, &v->value.z, data);
    break;
  case LP_VALUE_RATIONAL:
    rational_construct_copy(&v->value.q, data);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_construct_copy(&v->value.dy_q, data);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_copy(&v->value.a, data);
    break;
  }
}

void lp_value_construct_zero(lp_value_t* v) {
  v->type = LP_VALUE_INTEGER;
  integer_construct(&v->value.z);
}

void lp_value_construct_int(lp_value_t* v, long x) {
  v->type = LP_VALUE_INTEGER;
  integer_construct_from_int(lp_Z, &v->value.z, x);
}

void lp_value_construct_none(lp_value_t* v) {
  lp_value_construct(v, LP_VALUE_NONE, 0);
}

lp_value_t value_none = {
    .type = LP_VALUE_NONE
};

const lp_value_t* lp_value_none(void) {
  return &value_none;
}

lp_value_t value_minus_inf = {
    .type = LP_VALUE_MINUS_INFINITY
};

const lp_value_t* lp_value_minus_infinity(void) {
  return &value_minus_inf;
}

lp_value_t value_plus_inf = {
    .type = LP_VALUE_PLUS_INFINITY
};

const lp_value_t* lp_value_plus_infinity(void) {
  return &value_plus_inf;
}

void lp_value_assign_raw(lp_value_t* v, lp_value_type_t type, const void* data) {
  lp_value_destruct(v);
  lp_value_construct(v, type, data);
}

void lp_value_assign(lp_value_t* v, const lp_value_t* from) {
  if (v != from) {
    lp_value_destruct(v);
    lp_value_construct_copy(v, from);
  }
}

void lp_value_assign_zero(lp_value_t* v) {
  lp_value_destruct(v);
  lp_value_construct_zero(v);
}

void lp_value_swap(lp_value_t* v1, lp_value_t* v2) {
  lp_value_t tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

lp_value_t* lp_value_new(lp_value_type_t type, const void* data) {
  lp_value_t* result = malloc(sizeof(lp_value_t));
  lp_value_construct(result, type, data);
  return result;
}

lp_value_t* lp_value_new_copy(const lp_value_t* from) {
  lp_value_t* result = malloc(sizeof(lp_value_t));
  lp_value_construct_copy(result, from);
  return result;
}

void lp_value_construct_copy(lp_value_t* v, const lp_value_t* from) {
  switch(from->type) {
  case LP_VALUE_NONE:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    lp_value_construct(v, from->type, 0);
    break;
  case LP_VALUE_INTEGER:
    lp_value_construct(v, LP_VALUE_INTEGER, &from->value.z);
    break;
  case LP_VALUE_RATIONAL:
    lp_value_construct(v, LP_VALUE_RATIONAL, &from->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_value_construct(v, LP_VALUE_DYADIC_RATIONAL, &from->value.dy_q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_value_construct(v, LP_VALUE_ALGEBRAIC, &from->value.a);
    break;
  }
}

void lp_value_destruct(lp_value_t* v) {
  switch(v->type) {
  case LP_VALUE_NONE:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    break;
  case LP_VALUE_INTEGER:
    integer_destruct(&v->value.z);
    break;
  case LP_VALUE_RATIONAL:
    rational_destruct(&v->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_destruct(&v->value.dy_q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_destruct(&v->value.a);
    break;
  }
}

void lp_value_delete(lp_value_t* v) {
  lp_value_destruct(v);
  free(v);
}

// 1/2^20
#define LP_VALUE_APPROX_MIN_MAGNITUDE -20

void lp_value_approx(const lp_value_t* v, lp_rational_interval_t* out) {
  int size;

  lp_rational_interval_t approx;

  switch (v->type) {
  case LP_VALUE_NONE:
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    assert(0);
    break;
  case LP_VALUE_INTEGER: {
    lp_rational_t point;
    rational_construct_from_integer(&point, &v->value.z);
    lp_rational_interval_construct_point(&approx, &point);
    rational_destruct(&point);
    break;
  }
  case LP_VALUE_RATIONAL:
    lp_rational_interval_construct_point(&approx, &v->value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_rational_interval_construct_from_dyadic(&approx, &v->value.dy_q, 0, &v->value.dy_q, 0);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (lp_value_is_rational(v)) {
      lp_rational_t v_rat;
      rational_construct(&v_rat);
      lp_value_get_rational(v, &v_rat);
      lp_rational_interval_construct_point(&approx, &v_rat);
      rational_destruct(&v_rat);
    } else {
      // Make sure we're below the given size
      size = lp_dyadic_interval_size(&v->value.a.I);
      while (size > LP_VALUE_APPROX_MIN_MAGNITUDE) {
        size --;
        lp_algebraic_number_refine_const(&v->value.a);
      }
      lp_rational_interval_construct_from_dyadic_interval(&approx, &v->value.a.I);
    }
    break;
  }

  lp_rational_interval_swap(&approx, out);
  lp_rational_interval_destruct(&approx);
}

int lp_value_print(const lp_value_t* v, FILE* out) {
  int ret = 0;
  switch (v->type) {
  case LP_VALUE_NONE:
    ret += fprintf(out, "<null>");
    break;
  case LP_VALUE_PLUS_INFINITY:
    ret += fprintf(out, "+inf");
    break;
  case LP_VALUE_MINUS_INFINITY:
    ret += fprintf(out, "-inf");
    break;
  case LP_VALUE_INTEGER:
    ret += integer_print(&v->value.z, out);
    break;
  case LP_VALUE_RATIONAL:
    ret += rational_print(&v->value.q, out);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    ret += dyadic_rational_print(&v->value.dy_q, out);
    break;
  case LP_VALUE_ALGEBRAIC:
    ret += lp_algebraic_number_print(&v->value.a, out);
    break;
  }
  return ret;
}

int lp_value_sgn(const lp_value_t* v) {
  switch (v->type) {
  case LP_VALUE_NONE:
    assert(0);
    return 0;
  case LP_VALUE_PLUS_INFINITY:
    return 1;
  case LP_VALUE_MINUS_INFINITY:
    return -1;
  case LP_VALUE_INTEGER:
    return lp_integer_sgn(lp_Z, &v->value.z);
  case LP_VALUE_RATIONAL:
    return lp_rational_sgn(&v->value.q);
  case LP_VALUE_DYADIC_RATIONAL:
    return lp_dyadic_rational_sgn(&v->value.dy_q);
  case LP_VALUE_ALGEBRAIC:
    return lp_algebraic_number_sgn(&v->value.a);
  }
  return 0;
}

int lp_value_cmp(const lp_value_t* v1, const lp_value_t* v2) {

  if (trace_is_enabled("value::cmp")) {
    tracef("lp_value_cmp()\n")
    tracef("v1 = "); lp_value_print(v1, trace_out); tracef("\n");
    tracef("v2 = "); lp_value_print(v2, trace_out); tracef("\n");
  }

  if (v1 == v2) {
    return 0;
  }

  if (v1->type == v2->type) {
    lp_value_type_t type = v1->type;
    switch (type) {
    case LP_VALUE_NONE:
    case LP_VALUE_PLUS_INFINITY:
    case LP_VALUE_MINUS_INFINITY:
      return 0;
    case LP_VALUE_INTEGER:
      return lp_integer_cmp(lp_Z, &v1->value.z, &v2->value.z);
    case LP_VALUE_RATIONAL:
      return rational_cmp(&v1->value.q, &v2->value.q);
    case LP_VALUE_DYADIC_RATIONAL:
      return dyadic_rational_cmp(&v1->value.dy_q, &v2->value.dy_q);
    case LP_VALUE_ALGEBRAIC:
      return lp_algebraic_number_cmp(&v1->value.a, &v2->value.a);
    }
  }

  // Different types
  assert(v1->type != LP_VALUE_NONE && v2->type != LP_VALUE_NONE);
  if (v1->type == LP_VALUE_MINUS_INFINITY) {
    return -1;
  }
  if (v2->type == LP_VALUE_MINUS_INFINITY) {
    return 1;
  }
  if (v1->type == LP_VALUE_PLUS_INFINITY) {
    return 1;
  }
  if (v2->type == LP_VALUE_PLUS_INFINITY) {
    return -1;
  }

  // Make sure that the first one is bigger in the order int < dy_rat < rat < algebraic
  if (v1->type < v2->type) {
    return -lp_value_cmp(v2, v1);
  }

  switch (v1->type) {
  case LP_VALUE_DYADIC_RATIONAL:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      return dyadic_rational_cmp_integer(&v1->value.dy_q, &v2->value.z);
    default:
      assert(0);
    }
    break;
  case LP_VALUE_RATIONAL:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      return rational_cmp_integer(&v1->value.q, &v2->value.z);
    case LP_VALUE_DYADIC_RATIONAL:
      return rational_cmp_dyadic_rational(&v1->value.q, &v2->value.dy_q);
    default:
      assert(0);
    }
    break;
  case LP_VALUE_ALGEBRAIC:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      return lp_algebraic_number_cmp_integer(&v1->value.a, &v2->value.z);
    case LP_VALUE_DYADIC_RATIONAL:
      return lp_algebraic_number_cmp_dyadic_rational(&v1->value.a, &v2->value.dy_q);
    case LP_VALUE_RATIONAL:
      return lp_algebraic_number_cmp_rational(&v1->value.a, &v2->value.q);
    default:
      assert(0);
    }
    break;
  default:
    assert(0);
  }

  assert(0);

  return v1 == v2;
}

int lp_value_cmp_void(const void* v1, const void* v2) {
  return lp_value_cmp(v1, v2);
}

int lp_value_cmp_rational(const lp_value_t* v, const lp_rational_t* q) {

  int cmp  = 0;

  switch (v->type) {
  case LP_VALUE_PLUS_INFINITY:
    cmp = 1;
    break;
  case LP_VALUE_MINUS_INFINITY:
    cmp = -1;
    break;
  case LP_VALUE_INTEGER:
    cmp = -lp_rational_cmp_integer(q, &v->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    cmp = -lp_rational_cmp_dyadic_rational(q, &v->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    cmp = lp_rational_cmp(&v->value.q, q);
    break;
  case LP_VALUE_ALGEBRAIC:
    cmp = lp_algebraic_number_cmp_rational(&v->value.a, q);
    break;
  default:
    assert(0);
  }

  return cmp;
}

int lp_value_is_rational(const lp_value_t* v) {
  switch (v->type) {
  case LP_VALUE_INTEGER:
  case LP_VALUE_DYADIC_RATIONAL:
  case LP_VALUE_RATIONAL:
    return 1;
  case LP_VALUE_ALGEBRAIC:
    return lp_algebraic_number_is_rational(&v->value.a);
  default:
    return 0;
  }
}

int lp_value_is_integer(const lp_value_t* v) {
  switch (v->type) {
  case LP_VALUE_INTEGER:
    return 1;
  case LP_VALUE_DYADIC_RATIONAL:
    return lp_dyadic_rational_is_integer(&v->value.dy_q);
  case LP_VALUE_RATIONAL:
    return lp_rational_is_integer(&v->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    return lp_algebraic_number_is_integer(&v->value.a);
  default:
    return 0;
  }
}

int lp_value_is_infinity(const lp_value_t* v) {
  switch (v->type) {
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    return 1;
  default:
    return 0;
  }
}

void lp_value_ceiling(const lp_value_t* v, lp_integer_t* v_ceiling) {
  switch (v->type) {
  case LP_VALUE_INTEGER:
    lp_integer_assign(lp_Z, v_ceiling, &v->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_ceiling(&v->value.dy_q, v_ceiling);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_ceiling(&v->value.q, v_ceiling);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_ceiling(&v->value.a, v_ceiling);
    break;
  default:
    assert(0);
  }
}

void lp_value_floor(const lp_value_t* v, lp_integer_t* v_floor) {
  switch (v->type) {
  case LP_VALUE_INTEGER:
    lp_integer_assign(lp_Z, v_floor, &v->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_floor(&v->value.dy_q, v_floor);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_floor(&v->value.q, v_floor);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_floor(&v->value.a, v_floor);
    break;
  default:
    assert(0);
  }
}

void lp_value_get_rational(const lp_value_t* v, lp_rational_t* q) {
  lp_rational_t result;

  switch (v->type) {
  case LP_VALUE_INTEGER:
    rational_construct_from_integer(&result, &v->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    rational_construct_from_dyadic(&result, &v->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    rational_assign(q, &v->value.q);
    return;
  case LP_VALUE_ALGEBRAIC:
    if (lp_dyadic_interval_is_point(&v->value.a.I)) {
      // It's a point value, so we just get it
      lp_rational_construct_from_dyadic(&result, lp_dyadic_interval_get_point(&v->value.a.I));
    } else {
      const lp_upolynomial_t* v_poly = v->value.a.f;
      if (lp_upolynomial_degree(v_poly) == 1) {
        // p = ax + b = 0 => x = -b/a
        const lp_integer_t* b = lp_upolynomial_const_term(v_poly);
        const lp_integer_t* a = lp_upolynomial_lead_coeff(v_poly);
        if (b) {
          rational_construct_from_div(&result, b, a);
          rational_neg(&result, &result);
        } else {
          rational_construct(&result);
        }
      } else {
        assert(0);
      }
    }
    break;
  default:
    assert(0);
  }

  rational_swap(&result, q);
  rational_destruct(&result);
}

void lp_value_get_num(const lp_value_t* v, lp_integer_t* num) {
  assert(lp_value_is_rational(v));
  switch (v->type) {
  case LP_VALUE_INTEGER:
    integer_assign(lp_Z, num, &v->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_get_num(&v->value.dy_q, num);
    break;
  case LP_VALUE_RATIONAL:
    rational_get_num(&v->value.q, num);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (lp_dyadic_interval_is_point(&v->value.a.I)) {
      // It's a point value, so we just get it
      dyadic_rational_get_num(lp_dyadic_interval_get_point(&v->value.a.I), num);
    } else {
      const lp_upolynomial_t* v_poly = v->value.a.f;
      if (lp_upolynomial_degree(v_poly) == 1) {
        // p = ax + b = 0 => x = -b/a
        lp_rational_t value;
        const lp_integer_t* b = lp_upolynomial_const_term(v_poly);
        const lp_integer_t* a = lp_upolynomial_lead_coeff(v_poly);
        if (b) {
          rational_construct_from_div(&value, b, a);
          rational_neg(&value, &value);
        } else {
          rational_construct(&value);
        }
        rational_get_num(&value, num);
        rational_destruct(&value);
      } else {
        assert(0);
      }
    }
    break;
  default:
    assert(0);
  }
}

void lp_value_get_den(const lp_value_t* v, lp_integer_t* den) {
  assert(lp_value_is_rational(v));
  switch (v->type) {
  case LP_VALUE_INTEGER:
    integer_assign_int(lp_Z, den, 1);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    dyadic_rational_get_den(&v->value.dy_q, den);
    break;
  case LP_VALUE_RATIONAL:
    rational_get_den(&v->value.q, den);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (lp_dyadic_interval_is_point(&v->value.a.I)) {
      // It's a point value, so we just get it
      dyadic_rational_get_den(lp_dyadic_interval_get_point(&v->value.a.I), den);
    } else {
      const lp_upolynomial_t* v_poly = v->value.a.f;
      if (lp_upolynomial_degree(v_poly) == 1) {
        // p = ax + b = 0 => x = -b/a
        lp_rational_t value;
        const lp_integer_t* b = lp_upolynomial_const_term(v_poly);
        const lp_integer_t* a = lp_upolynomial_lead_coeff(v_poly);
        if (b) {
          rational_construct_from_div(&value, b, a);
          rational_neg(&value, &value);
        } else {
          rational_construct(&value);
        }
        rational_get_den(&value, den);
        rational_destruct(&value);
      } else {
        assert(0);
      }
    }
    break;
  default:
    assert(0);
  }
}

char* lp_value_to_string(const lp_value_t* v) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_value_print(v, f);
  fclose(f);
  return str;
}

void lp_value_get_value_between(const lp_value_t* a, int a_strict, const lp_value_t* b, int b_strict, lp_value_t* v) {

  if (trace_is_enabled("value::get_value_between")) {
    tracef("lp_value_get_value_between()\n")
    tracef("a = "); lp_value_print(a, trace_out); tracef(", a_strict = %d\n", a_strict);
    tracef("b = "); lp_value_print(b, trace_out); tracef(", b_strict = %d\n", b_strict);
  }

  // Compare a and b. If they are equal the only option is the value they hold.
  // If a > b then we swap them
  int cmp = lp_value_cmp(a, b);
  if (cmp == 0) {
    // Same, we're done
    assert(!a_strict && !b_strict);
    lp_value_assign(v, a);
    return;
  } else if (cmp > 0) {
    const lp_value_t* tmp = a; a = b; b = tmp;
    int strict_tmp = a_strict; a_strict = b_strict; b_strict = strict_tmp;
  }

  // If the whole R we're done, just pick 0
  if (a->type == LP_VALUE_MINUS_INFINITY && b->type == LP_VALUE_PLUS_INFINITY) {
     // (-inf, +inf), just pick 0
     lp_integer_t zero;
     lp_integer_construct(&zero);
     lp_value_assign_raw(v, LP_VALUE_INTEGER, &zero);
     lp_integer_destruct(&zero);
     return;
   }

  // We have a < b, at least one of them is not infinity, and comparison ensures
  // that the algebraic intervals will be disjoint

  // Get rational values a_ub and b_lb such that a <= a_ub <= b_lb <= b
  lp_rational_t a_ub, b_lb;
  int a_inf = 0, b_inf = 0;
  int a_ub_strict = a_strict;
  int b_lb_strict = b_strict;


  switch (a->type) {
  case LP_VALUE_MINUS_INFINITY:
    a_inf = 1;
    break;
  case LP_VALUE_INTEGER:
    rational_construct_from_integer(&a_ub, &a->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    rational_construct_from_dyadic(&a_ub, &a->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    rational_construct_copy(&a_ub, &a->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (lp_value_is_rational(a)) {
      rational_construct(&a_ub);
      lp_value_get_rational(a, &a_ub);
    } else {
      // Get the upper bound of the interval as a_ub
      rational_construct_from_dyadic(&a_ub, &a->value.a.I.b);
      // Algebaic bound is strict so a_ub can be picked
      a_ub_strict = 0;
    }
    break;
  default:
    assert(0);
  }

  switch (b->type) {
  case LP_VALUE_PLUS_INFINITY:
    b_inf = 1;
    break;
  case LP_VALUE_INTEGER:
    rational_construct_from_integer(&b_lb, &b->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    rational_construct_from_dyadic(&b_lb, &b->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    rational_construct_copy(&b_lb, &b->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (lp_value_is_rational(b)) {
      rational_construct(&b_lb);
      lp_value_get_rational(b, &b_lb);
    } else {
      // Get the lower bound of the interval as b_lb
      rational_construct_from_dyadic(&b_lb, &b->value.a.I.a);
      // Algebraic bound is strict so b_lb can be picked
      b_lb_strict = 0;
    }
    break;
  default:
    assert(0);
  }

  assert(!a_inf || !b_inf);

  // We have rational values a_ub and b_ub such that a <= a_ub <= b_lb <= b
  // or one of them is infinity

  if (a_inf) {
    // -inf < b_lb <= b
    // just pick the value to be floor(b_lb)-1
    lp_integer_t result;
    lp_integer_construct(&result);
    rational_floor(&b_lb, &result);
    lp_integer_dec(lp_Z, &result);
    lp_value_assign_raw(v, LP_VALUE_INTEGER, &result);
    lp_integer_destruct(&result);
    lp_rational_destruct(&b_lb);
    return;
  }

  if (b_inf) {
    // a <= a_ub < +inf
    // just pick the value to be ceil(a_ub)+1
    lp_integer_t result;
    lp_integer_construct(&result);
    rational_ceiling(&a_ub, &result);
    lp_integer_inc(lp_Z, &result);
    lp_value_assign_raw(v, LP_VALUE_INTEGER, &result);
    lp_integer_destruct(&result);
    lp_rational_destruct(&a_ub);
    return;
  }

  // We have rational values a_ub and b_ub such that a <= a_ub <= b_lb <= b
  // Both a_ub and b_lb are constructed

  // If a_ub == b_lb, this is due to algebraic number intervals, so refine once more
  cmp = rational_cmp(&a_ub, &b_lb);
  if (cmp == 0) {
    assert(!lp_value_is_rational(a) || !lp_value_is_rational(b));
    if (!lp_value_is_rational(a)) {
      assert(a->type == LP_VALUE_ALGEBRAIC);
      lp_algebraic_number_refine_const(&a->value.a);
    }
    if (!lp_value_is_rational(b)) {
      assert(b->type == LP_VALUE_ALGEBRAIC);
      lp_algebraic_number_refine_const(&b->value.a);
    }
    lp_value_get_value_between(a, a_strict, b, b_strict, v);
  } else {

    // To be constructed to the value
    lp_rational_t result;

    // Get the smallest integer interval around [a_ub, b_lb] and refine

    // a <= a_ub < b_lb <= b
    // m = (a_ub + b_lb)/2 so a < m < b
    lp_rational_t m;
    rational_construct(&m);
    rational_add(&m, &a_ub, &b_lb);
    rational_div_2exp(&m, &m, 1);

    // floor(m) <= m <= ceil(m)
    lp_integer_t m_floor, m_ceil;
    integer_construct(&m_floor);
    rational_floor(&m, &m_floor);
    integer_construct_copy(lp_Z, &m_ceil, &m_floor);
    integer_inc(lp_Z, &m_ceil);

    if (trace_is_enabled("value::get_value_between")) {
      tracef("a_ub = ");
      lp_rational_print(&a_ub, trace_out);
      tracef("\n");
      tracef("b_ub = ");
      lp_rational_print(&b_lb, trace_out);
      tracef("\n");
      tracef("m = ");
      lp_rational_print(&m, trace_out);
      tracef("\n");
      tracef("m_floor = ");
      lp_integer_print(&m_floor, trace_out);
      tracef("\n");
      tracef("m_ceil = ");
      lp_integer_print(&m_ceil, trace_out);
      tracef("\n");
    }

    // If a_ub < m_floor (or equal and bound allows it), we can take this value
    cmp = lp_rational_cmp_integer(&a_ub, &m_floor);
    if (cmp < 0 || (cmp == 0 && !a_ub_strict)) {
      lp_rational_construct_from_integer(&result, &m_floor);
    } else {
      // If m_ceil < b_lb (or equal and bound allows it), we can take this value
      cmp = lp_rational_cmp_integer(&b_lb, &m_ceil);
      if (cmp > 0 || (cmp == 0 && !b_lb_strict)) {
        lp_rational_construct_from_integer(&result, &m_ceil);
      } else {

        lp_rational_t lb, ub;
        rational_construct_from_integer(&lb, &m_floor);
        rational_construct_from_integer(&ub, &m_ceil);

        // We have to do the search
        for (;;) {

          // always: lb < a_ub <= b_lb < ub
          rational_add(&m, &lb, &ub);
          rational_div_2exp(&m, &m, 1);

          // lb < m < a_ub => move lb to m
          cmp = rational_cmp(&a_ub, &m);
          // if ((a_strict && cmp >= 0) || (!a_strict && cmp > 0)) {
          if (cmp >= 0) {
            rational_swap(&m, &lb);
            continue;
          }

          // b_lb < m < ub => move ub to m
          cmp = rational_cmp(&m, &b_lb);
          // if ((b_strict && cmp >= 0) || (!b_strict && cmp > 0)) {
          if (cmp >= 0) {
            rational_swap(&ub, &m);
            continue;
          }

          // Got it l <= m <= u
          rational_construct_copy(&result, &m);
          break;
        }

        rational_destruct(&lb);
        rational_destruct(&ub);
      }
    }

    lp_value_assign_raw(v, LP_VALUE_RATIONAL, &result);

    rational_destruct(&result);
    integer_destruct(&m_ceil);
    integer_destruct(&m_floor);
    rational_destruct(&m);
  }

  rational_destruct(&a_ub);
  rational_destruct(&b_lb);

  if (trace_is_enabled("value::get_value_between")) {
    tracef("lp_value_get_value_between() => "); lp_value_print(v, trace_out); tracef("\n");
  }

}

int lp_value_get_distance_size_approx(const lp_value_t* lower, const lp_value_t* upper) {

  assert(lp_value_cmp(lower, upper) < 0);

  if (lower->type == LP_VALUE_MINUS_INFINITY) {
    return INT_MAX;
  }

  if (upper->type == LP_VALUE_PLUS_INFINITY) {
    return INT_MAX;
  }

  lp_rational_t lower_approx, upper_approx;
  rational_construct(&lower_approx);
  rational_construct(&upper_approx);

  // Get lower bound approximation
  if (lp_value_is_rational(lower)) {
    lp_value_get_rational(lower, &lower_approx);
  } else {
    assert(lower->type == LP_VALUE_ALGEBRAIC);
    assert(!lower->value.a.I.is_point);
    lp_algebraic_number_get_rational_midpoint(&lower->value.a, &lower_approx);
  }

  if (lp_value_is_rational(upper)) {
    lp_value_get_rational(upper, &upper_approx);
  } else {
    assert(upper->type == LP_VALUE_ALGEBRAIC);
    assert(!upper->value.a.I.is_point);
    lp_algebraic_number_get_rational_midpoint(&upper->value.a, &upper_approx);
  }

  // Get the distance
  lp_rational_t* m = &lower_approx;
  lp_rational_sub(m, &upper_approx, &lower_approx);

  // The denominator and numerator
  lp_integer_t num, den;
  integer_construct(&num);
  integer_construct(&den);
  rational_get_num(m, &num);
  rational_get_den(m, &den);

  // Size = log(num/den) = log(num) - log(den)
  int size = ((int)integer_log2_abs(&num)) - ((int)integer_log2_abs(&den)) + 1;

  integer_destruct(&num);
  integer_destruct(&den);
  rational_destruct(&lower_approx);
  rational_destruct(&upper_approx);

  return size;
}

size_t lp_value_hash(const lp_value_t* v) {
  return lp_value_hash_approx(v, 0);
}

size_t lp_value_hash_approx(const lp_value_t* v, unsigned precision) {
  switch (v->type) {
  case LP_VALUE_NONE:
    return 0;
  case LP_VALUE_PLUS_INFINITY:
    return SIZE_MAX-1;
  case LP_VALUE_MINUS_INFINITY:
    return SIZE_MAX;
  case LP_VALUE_INTEGER:
    return lp_integer_hash(&v->value.z);
  case LP_VALUE_DYADIC_RATIONAL:
    return lp_dyadic_rational_hash_approx(&v->value.dy_q, precision);
  case LP_VALUE_RATIONAL:
    return lp_rational_hash_approx(&v->value.q, precision);
  case LP_VALUE_ALGEBRAIC:
    return lp_algebraic_number_hash_approx(&v->value.a, precision);
  default:
    assert(0);
    return 0;
  }
}


double lp_value_to_double(const lp_value_t* v) {
  switch (v->type) {
  case LP_VALUE_NONE:
    return 0;
  case LP_VALUE_PLUS_INFINITY:
    return INFINITY;
  case LP_VALUE_MINUS_INFINITY:
    return -INFINITY;
  case LP_VALUE_INTEGER:
    return lp_integer_to_double(&v->value.z);
  case LP_VALUE_DYADIC_RATIONAL:
    return lp_dyadic_rational_to_double(&v->value.dy_q);
  case LP_VALUE_RATIONAL:
    return lp_rational_to_double(&v->value.q);
  case LP_VALUE_ALGEBRAIC:
    return lp_algebraic_number_to_double(&v->value.a);
  default:
    assert(0);
    return 0;
  }
}

/**
 * Case v1 and v2 to the same type.
 *
 * @param v1, v2 the values to case
 * @param v1_tmp the value to use for allocating a new value
 * @param v1_to_use, v2_to_use the cast values to use
 */
int lp_value_to_same_type(const lp_value_t* v1, const lp_value_t* v2,
    lp_value_t* v1_new, lp_value_t* v2_new,
    const lp_value_t** v1_to_use, const lp_value_t** v2_to_use) {

  // trivial
  if (v1->type == v2->type) {
    *v1_to_use = v1;
    *v2_to_use = v2;
    return 1;
  }

  lp_rational_t tmp_rat;
  lp_dyadic_rational_t tmp_dy;
  lp_algebraic_number_t tmp_a;

  switch (v1->type) {
  case LP_VALUE_INTEGER:
    switch (v2->type) {
    case LP_VALUE_DYADIC_RATIONAL:
      // v1: integer
      // v2: dyadic rational
      lp_dyadic_rational_construct_from_integer(&tmp_dy, &v1->value.z);
      lp_value_construct(v1_new, LP_VALUE_DYADIC_RATIONAL, &tmp_dy);
      lp_dyadic_rational_destruct(&tmp_dy);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    case LP_VALUE_RATIONAL:
      // v1: integer
      // v2: rational
      lp_rational_construct_from_integer(&tmp_rat, &v1->value.z);
      lp_value_construct(v1_new, LP_VALUE_RATIONAL, &tmp_rat);
      lp_rational_destruct(&tmp_rat);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    case LP_VALUE_ALGEBRAIC:
      // v1: integer
      // v2: algebraic
      lp_algebraic_number_construct_from_integer(&tmp_a, &v1->value.z);
      lp_value_construct(v1_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    default:
      // unsupported
      return 0;
    }
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      // v1: dyadic rational
      // v2: integer
      lp_dyadic_rational_construct_from_integer(&tmp_dy, &v2->value.z);
      lp_value_construct(v2_new, LP_VALUE_DYADIC_RATIONAL, &tmp_dy);
      lp_dyadic_rational_destruct(&tmp_dy);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    case LP_VALUE_RATIONAL:
      // v1: dyadic rational
      // v2: rational
      lp_rational_construct_from_dyadic(&tmp_rat, &v1->value.dy_q);
      lp_value_construct(v1_new, LP_VALUE_RATIONAL, &tmp_rat);
      lp_rational_destruct(&tmp_rat);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    case LP_VALUE_ALGEBRAIC:
      // v1: dyadic rational
      // v2: algebraic
      lp_algebraic_number_construct_from_dyadic_rational(&tmp_a, &v1->value.dy_q);
      lp_value_construct(v1_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    default:
      // unsupported
      return 0;
    }
    break;
  case LP_VALUE_RATIONAL:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      // v1: rational
      // v2: integer
      lp_rational_construct_from_integer(&tmp_rat, &v2->value.z);
      lp_value_construct(v2_new, LP_VALUE_RATIONAL, &tmp_rat);
      lp_rational_destruct(&tmp_rat);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    case LP_VALUE_DYADIC_RATIONAL:
      // v1: rational
      // v2: dyadic rational
      lp_rational_construct_from_dyadic(&tmp_rat, &v2->value.dy_q);
      lp_value_construct(v2_new, LP_VALUE_RATIONAL, &tmp_rat);
      lp_rational_destruct(&tmp_rat);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    case LP_VALUE_ALGEBRAIC:
      // v1: rational
      // v2: algebraic_number
      lp_algebraic_number_construct_from_rational(&tmp_a, &v1->value.q);
      lp_value_construct(v2_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1_new;
      *v2_to_use = v2;
      break;
    default:
      // unsupported
      return 0;
    }
    break;
  case LP_VALUE_ALGEBRAIC:
    switch (v2->type) {
    case LP_VALUE_INTEGER:
      // v1: algebraic number
      // v2: integer
      lp_algebraic_number_construct_from_integer(&tmp_a, &v2->value.z);
      lp_value_construct(v2_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    case LP_VALUE_DYADIC_RATIONAL:
      // v1: algebraic number
      // v2: dyadic rational
      lp_algebraic_number_construct_from_dyadic_rational(&tmp_a, &v2->value.dy_q);
      lp_value_construct(v2_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    case LP_VALUE_RATIONAL:
      // v1: algebraic number
      // v2: ratioal number
      lp_algebraic_number_construct_from_rational(&tmp_a, &v2->value.q);
      lp_value_construct(v2_new, LP_VALUE_ALGEBRAIC, &tmp_a);
      lp_algebraic_number_destruct(&tmp_a);
      *v1_to_use = v1;
      *v2_to_use = v2_new;
      break;
    default:
      // unsupported
      return 0;
    }

  default:
    // unsupported
    return 0;
  }

  return 1;
}

void lp_value_add(lp_value_t* sum, const lp_value_t* a, const lp_value_t* b) {

  lp_value_t a_new, b_new;
  const lp_value_t *a_to_use;
  const lp_value_t *b_to_use;

  // Check for infinities
  if (a->type == LP_VALUE_PLUS_INFINITY) {
    if (b->type != LP_VALUE_MINUS_INFINITY) {
      lp_value_assign_raw(sum, LP_VALUE_PLUS_INFINITY, 0);
    } else {
      assert(0);
    }
    return;
  } else if (a->type == LP_VALUE_MINUS_INFINITY) {
    if (b->type != LP_VALUE_PLUS_INFINITY) {
      lp_value_assign_raw(sum, LP_VALUE_MINUS_INFINITY, 0);
    } else {
      assert(0);
    }
    return;
  } else if (b->type == LP_VALUE_PLUS_INFINITY) {
    if (a->type != LP_VALUE_MINUS_INFINITY) {
      lp_value_assign_raw(sum, LP_VALUE_PLUS_INFINITY, 0);
    } else {
      assert(0);
    }
    return;
  } else if (b->type == LP_VALUE_MINUS_INFINITY) {
    if (a->type != LP_VALUE_PLUS_INFINITY) {
      lp_value_assign_raw(sum, LP_VALUE_MINUS_INFINITY, 0);
    } else {
      assert(0);
    }
    return;
  }

  // All other cases, cast and run the operation
  int ret = lp_value_to_same_type(a, b, &a_new, &b_new, &a_to_use, &b_to_use);
  (void) ret;
  assert(ret);

  lp_value_t result;

  result.type = a_to_use->type;
  switch (result.type) {
  case LP_VALUE_INTEGER:
    lp_integer_construct(&result.value.z);
    lp_integer_add(lp_Z, &result.value.z, &a_to_use->value.z, &b_to_use->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_construct(&result.value.dy_q);
    lp_dyadic_rational_add(&result.value.dy_q, &a_to_use->value.dy_q, &b_to_use->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&result.value.q);
    lp_rational_add(&result.value.q, &a_to_use->value.q, &b_to_use->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_zero(&result.value.a);
    lp_algebraic_number_add(&result.value.a, &a_to_use->value.a, &b_to_use->value.a);
    break;
  default:
    assert(0);
    break;
  }

  if (a_to_use != a) {
    lp_value_destruct(a_to_use);
  }
  if (b_to_use != b) {
    lp_value_destruct(b_to_use);
  }

  lp_value_swap(sum, &result);
  lp_value_destruct(&result);
}

void lp_value_sub(lp_value_t* sub, const lp_value_t* a, const lp_value_t* b) {
  lp_value_t b_neg;
  lp_value_construct_none(&b_neg);
  lp_value_neg(&b_neg, b);
  lp_value_add(sub, a, &b_neg);
  lp_value_destruct(&b_neg);
}

void lp_value_neg(lp_value_t* neg, const lp_value_t* a) {

  lp_value_t result;

  result.type = a->type;
  switch(a->type) {
  case LP_VALUE_NONE:
    break;
  case LP_VALUE_INTEGER:
    lp_integer_construct(&result.value.z);
    lp_integer_neg(lp_Z, &result.value.z, &a->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_construct(&result.value.dy_q);
    lp_dyadic_rational_neg(&result.value.dy_q, &a->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&result.value.q);
    lp_rational_neg(&result.value.q, &a->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_zero(&result.value.a);
    lp_algebraic_number_neg(&result.value.a, &a->value.a);
    break;
  case LP_VALUE_PLUS_INFINITY:
    result.type = LP_VALUE_MINUS_INFINITY;
    break;
  case LP_VALUE_MINUS_INFINITY:
    result.type = LP_VALUE_PLUS_INFINITY;
    break;
  }

  lp_value_swap(neg, &result);
  lp_value_destruct(&result);
}

void lp_value_mul(lp_value_t* mul, const lp_value_t* a, const lp_value_t* b) {

  lp_value_t a_new, b_new;
  const lp_value_t *a_to_use;
  const lp_value_t *b_to_use;

  // Check for infinities
  if (lp_value_is_infinity(a) || lp_value_is_infinity(b)) {
    int sgn_a = lp_value_sgn(a);
    int sgn_b = lp_value_sgn(b);
    int sgn_mul = sgn_a * sgn_b;
    if (sgn_mul > 0) {
      lp_value_assign_raw(mul, LP_VALUE_PLUS_INFINITY, 0);
    } else if (sgn_mul) {
      lp_value_assign_raw(mul, LP_VALUE_MINUS_INFINITY, 0);
    } else {
      assert(0);
    }
    return;
  }

  // All other cases, cast and run the operation
  int ret = lp_value_to_same_type(a, b, &a_new, &b_new, &a_to_use, &b_to_use);
  (void) ret;
  assert(ret);

  lp_value_t result;

  result.type = a_to_use->type;
  switch (result.type) {
  case LP_VALUE_INTEGER:
    lp_integer_construct(&result.value.z);
    lp_integer_mul(lp_Z, &result.value.z, &a_to_use->value.z, &b_to_use->value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_construct(&result.value.dy_q);
    lp_dyadic_rational_mul(&result.value.dy_q, &a_to_use->value.dy_q, &b_to_use->value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&result.value.q);
    lp_rational_mul(&result.value.q, &a_to_use->value.q, &b_to_use->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_zero(&result.value.a);
    lp_algebraic_number_mul(&result.value.a, &a_to_use->value.a, &b_to_use->value.a);
    break;
  default:
    assert(0);
    break;
  }

  if (a_to_use != a) {
    lp_value_destruct(a_to_use);
  }
  if (b_to_use != b) {
    lp_value_destruct(b_to_use);
  }

  lp_value_swap(mul, &result);
  lp_value_destruct(&result);
}

void lp_value_inv(lp_value_t* inv, const lp_value_t* a) {

  lp_value_t result;

  switch(inv->type) {
  case LP_VALUE_NONE:
    break;
  case LP_VALUE_INTEGER:
    result.type = LP_VALUE_RATIONAL;
    lp_rational_construct_from_integer(&result.value.q, &a->value.z);
    lp_rational_inv(&result.value.q, &result.value.q);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    result.type = LP_VALUE_RATIONAL;
    lp_rational_construct_from_dyadic(&result.value.q, &a->value.dy_q);
    lp_rational_inv(&result.value.q, &result.value.q);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&result.value.q);
    lp_rational_inv(&result.value.q, &a->value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    result.type = LP_VALUE_RATIONAL;
    lp_algebraic_number_construct_zero(&result.value.a);
    lp_algebraic_number_inv(&result.value.a, &a->value.a);
    break;
  case LP_VALUE_PLUS_INFINITY:
  case LP_VALUE_MINUS_INFINITY:
    lp_value_construct_zero(&result);
    break;
  }

  lp_value_swap(inv, &result);
  lp_value_destruct(&result);
}

void lp_value_div(lp_value_t* div, const lp_value_t* a, const lp_value_t* b) {
  lp_value_t b_inv;
  lp_value_construct_none(&b_inv);
  lp_value_inv(&b_inv, b);
  lp_value_mul(div, a, &b_inv);
  lp_value_destruct(&b_inv);
}

void lp_value_pow(lp_value_t* pow, const lp_value_t* a, unsigned n) {

  lp_value_t result;

  result.type = a->type;
  switch(a->type) {
  case LP_VALUE_NONE:
    break;
  case LP_VALUE_INTEGER:
    lp_integer_construct(&result.value.z);
    lp_integer_pow(lp_Z, &result.value.z, &a->value.z, n);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_construct(&result.value.dy_q);
    lp_dyadic_rational_pow(&result.value.dy_q, &a->value.dy_q, n);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&result.value.q);
    lp_rational_pow(&result.value.q, &a->value.q, n);
    break;
  case LP_VALUE_ALGEBRAIC:
    lp_algebraic_number_construct_zero(&result.value.a);
    lp_algebraic_number_pow(&result.value.a, &a->value.a, n);
    break;
  case LP_VALUE_PLUS_INFINITY:
    result.type = LP_VALUE_MINUS_INFINITY;
    break;
  case LP_VALUE_MINUS_INFINITY:
    if (n % 2) {
      result.type = LP_VALUE_MINUS_INFINITY;
    } else {
      result.type = LP_VALUE_PLUS_INFINITY;
    }
    break;
  }

  lp_value_swap(pow, &result);
  lp_value_destruct(&result);
}

