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

#include "interval/arithmetic.h"

#include "number/rational.h"
#include "number/dyadic_rational.h"
#include "number/value.h"

int rational_interval_endpoint_lt(const lp_rational_t* a, int a_open, const lp_rational_t* b, int b_open) {
  int cmp = rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

int dyadic_interval_endpoint_lt(const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t*b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

int lp_interval_endpoint_lt(const lp_value_t* a, int a_open, const lp_value_t* b, int b_open) {
  int cmp = lp_value_cmp(a, b);
  if (cmp == 0) {
    return (!a_open && b_open);
  } else {
    return cmp < 0;
  }
}

void rational_interval_add(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {

  if (I1->is_point && I2->is_point) {
    if (!S->is_point) {
      rational_destruct(&S->b);
    }
    rational_add(&S->a, &I1->a, &I2->a);
    S->b_open = S->a_open = 0;
    S->is_point = 1;
    return;
  }

  if (I2->is_point) {
    // Reuse symmetry
    rational_interval_add(S, I2, I1);
    return;
  }

  lp_rational_interval_t result;

  if (I1->is_point) {
    // Just shift by I1->a (I2 is not a point)
    lp_rational_interval_construct_copy(&result, I2);
    rational_add(&result.a, &result.a, &I1->a);
    rational_add(&result.b, &result.b, &I1->a);
  } else {
    // [a, b] + [c, d] = [a + c, b + d]
    rational_construct(&result.a);
    rational_construct(&result.b);
    rational_add(&result.a, &I1->a, &I2->a);
    rational_add(&result.b, &I1->b, &I2->b);
    result.a_open = I1->a_open || I2->a_open;
    result.b_open = I1->b_open || I2->b_open;
    result.is_point = 0;
  }

  lp_rational_interval_swap(&result, S);
  lp_rational_interval_destruct(&result);
}

void dyadic_interval_add(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {

  if (I1->is_point && I2->is_point) {
    if (!S->is_point) {
      dyadic_rational_destruct(&S->b);
    }
    dyadic_rational_add(&S->a, &I1->a, &I2->a);
    S->b_open = S->a_open = 0;
    S->is_point = 1;
    return;
  }

  if (I2->is_point) {
    // Reuse symmetry
    dyadic_interval_add(S, I2, I1);
    return;
  }

  lp_dyadic_interval_t result;

  if (I1->is_point) {
    // Just shift by I1->a
    lp_dyadic_interval_construct_copy(&result, I2);
    dyadic_rational_add(&result.a, &result.a, &I1->a);
    dyadic_rational_add(&result.b, &result.b, &I1->a);
  } else {
    // Both non-points
    // [a, b] + [c, d] = [a + c, b + d]
    dyadic_rational_construct(&result.a);
    dyadic_rational_construct(&result.b);
    dyadic_rational_add(&result.a, &I1->a, &I2->a);
    dyadic_rational_add(&result.b, &I1->b, &I2->b);
    result.a_open = I1->a_open || I2->a_open;
    result.b_open = I1->b_open || I2->b_open;
    result.is_point = 0;
  }

  lp_dyadic_interval_swap(&result, S);
  lp_dyadic_interval_destruct(&result);
}

static
int lp_value_add_approx(const lp_value_t* v1, const lp_value_t* v2, lp_value_t* lb, lp_value_t* ub) {

  assert(v1 != lb && v1 != ub);
  assert(v2 != lb && v2 != ub);
  assert(lb != ub);

  lp_integer_t add_int;
  lp_rational_t add_rat;
  lp_dyadic_rational_t add_dy;
  lp_dyadic_interval_t add_dy_interval;

  int is_point = 1;

  assert(v1->type != LP_VALUE_NONE);
  assert(v2->type != LP_VALUE_NONE);

  if (v1->type == v2->type) {
    switch (v1->type) {
    case LP_VALUE_MINUS_INFINITY:
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_MINUS_INFINITY, 0);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_MINUS_INFINITY, 0);
      }
      break;
    case LP_VALUE_PLUS_INFINITY:
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_PLUS_INFINITY, 0);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_PLUS_INFINITY, 0);
      }
      break;
    case LP_VALUE_INTEGER:
      lp_integer_construct(&add_int);
      lp_integer_add(lp_Z, &add_int, &v1->value.z, &v2->value.z);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_INTEGER, &add_int);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_INTEGER, &add_int);
      }
      lp_integer_destruct(&add_int);
      break;
    case LP_VALUE_RATIONAL:
      lp_rational_construct(&add_rat);
      lp_rational_add(&add_rat, &v1->value.q, &v2->value.q);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_RATIONAL, &add_rat);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_RATIONAL, &add_rat);
      }
      lp_rational_destruct(&add_rat);
      break;
    case LP_VALUE_DYADIC_RATIONAL:
      lp_dyadic_rational_construct(&add_dy);
      lp_dyadic_rational_add(&add_dy, &v1->value.dy_q, &v2->value.dy_q);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &add_dy);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &add_dy);
      }
      lp_dyadic_rational_destruct(&add_dy);
      break;
    case LP_VALUE_ALGEBRAIC:
      if (v1->value.a.I.is_point && v2->value.a.I.is_point) {
        lp_dyadic_rational_construct(&add_dy);
        lp_dyadic_rational_add(&add_dy, &v1->value.a.I.a, &v2->value.a.I.a);
        if (lb) {
          lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &add_dy);
        }
        if (ub) {
          lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &add_dy);
        }
        lp_dyadic_rational_destruct(&add_dy);
      } else {
        lp_dyadic_interval_construct_zero(&add_dy_interval);
        dyadic_interval_add(&add_dy_interval, &v1->value.a.I, &v2->value.a.I);
        if (lb) {
          lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &add_dy_interval.a);
        }
        if (ub) {
          lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &add_dy_interval.b);
        }
        is_point = 0;
      }
      break;
    case LP_VALUE_NONE:
      assert(0);
      break;
    }
  } else {

    lp_value_t v1_tmp, v2_tmp;
    const lp_value_t* v1_to_use;
    const lp_value_t* v2_to_use;

    // Cast to the same type
    int cast = lp_value_to_same_type(v1, v2, &v1_tmp, &v2_tmp, &v1_to_use, &v2_to_use);
    if (cast) {
      // Cast successful, just add as usual
      lp_value_add_approx(v1_to_use, v2_to_use, lb, ub);
      // If a fesh value was used, delete it
      if (v1_to_use != v1) {
        lp_value_destruct(&v1_tmp);
      }
      if (v2_to_use != v2) {
        lp_value_destruct(&v2_tmp);
      }
    } else {
      // Cast failed, check the cases
      int undefined = 0;
      lp_value_type_t result = LP_VALUE_NONE;

      // Check infinity cases
      if (v1->type == LP_VALUE_MINUS_INFINITY) {
        if (v2->type != LP_VALUE_PLUS_INFINITY) {
          result = LP_VALUE_MINUS_INFINITY;
        } else {
          undefined = 1;
        }
      }
      if (v1->type == LP_VALUE_PLUS_INFINITY) {
        if (v2->type != LP_VALUE_MINUS_INFINITY) {
          result = LP_VALUE_PLUS_INFINITY;
        } else {
          undefined = 1;
        }
      }
      if (v2->type == LP_VALUE_MINUS_INFINITY) {
        if (v1->type != LP_VALUE_PLUS_INFINITY) {
          result = LP_VALUE_MINUS_INFINITY;
        } else {
          undefined = 1;
        }
      }
      if (v2->type == LP_VALUE_PLUS_INFINITY) {
        if (v1->type != LP_VALUE_MINUS_INFINITY) {
          result = LP_VALUE_PLUS_INFINITY;
        } else {
          undefined = 1;
        }
      }

      // If infinity handled, we're done
      if (result != LP_VALUE_NONE) {
        if (lb) {
          lp_value_assign_raw(lb, result, 0);
        }
        if (ub) {
          lp_value_assign_raw(ub, result, 0);
        }
        return 1;
      }
      // If undefined, return NONE
      if (undefined) {
        if (lb) {
          lp_value_assign_raw(lb, LP_VALUE_NONE, 0);
        }
        if (ub) {
          lp_value_assign_raw(ub, LP_VALUE_NONE, 0);
        }
        return 0;
      }

      // OK, now deal with algebraic numbers, assume v1 is algebraic
      if (v2->type == LP_VALUE_ALGEBRAIC) {
        const lp_value_t* tmp = v1; v1 = v2; v2 = tmp;
      }

      assert(v1->type == LP_VALUE_ALGEBRAIC);
      assert(v2->type != LP_VALUE_MINUS_INFINITY);
      assert(v2->type != LP_VALUE_PLUS_INFINITY);
      assert(v2->type != LP_VALUE_ALGEBRAIC);

      // v1 is of the form (dy1, dy2), v2 is a proper non-algebraic value
      // we add them into (dy1 + v2, dy2 + v2)

      lp_value_t v1_lb, v1_ub;
      const lp_dyadic_rational_t* a = &v1->value.a.I.a;
      const lp_dyadic_rational_t* b = v1->value.a.I.is_point ? a : &v1->value.a.I.b;
      lp_value_construct(&v1_lb, LP_VALUE_DYADIC_RATIONAL, a);
      lp_value_construct(&v1_ub, LP_VALUE_DYADIC_RATIONAL, b);
      lp_value_add_approx(&v1_lb, v2, lb, 0);
      lp_value_add_approx(&v1_ub, v2, ub, 0);
      is_point = 0;
      lp_value_destruct(&v1_lb);
      lp_value_destruct(&v1_ub);
    }
  }

  return is_point;
}

void lp_interval_add(lp_interval_t* add, const lp_interval_t* I1, const lp_interval_t* I2) {

  lp_interval_t result;
  lp_interval_construct_full(&result);

  if (I1->is_point && I2->is_point) {
    result.is_point = lp_value_add_approx(&I1->a, &I2->a, &result.a, &result.b);
    if (result.is_point) {
      lp_value_destruct(&result.b);
    }
    result.b_open = result.a_open = !result.is_point;
  } else {
    // [a, b] + [c, d] = [a + c, b + d]
    const lp_value_t* I1_a = &I1->a;
    const lp_value_t* I1_b = I1->is_point ? I1_a : &I1->b;
    const lp_value_t* I2_a = &I2->a;
    const lp_value_t* I2_b = I2->is_point ? I2_a : &I2->b;
    int a_point = lp_value_add_approx(I1_a, I2_a, &result.a, 0);
    int b_point = lp_value_add_approx(I1_b, I2_b, 0, &result.b);
    result.a_open = I1->a_open || I2->a_open || !a_point;
    result.b_open = I1->b_open || I2->b_open || !b_point;
    result.is_point = 0;
  }

  lp_interval_swap(add, &result);
  lp_interval_destruct(&result);
}

void rational_interval_neg(lp_rational_interval_t* N, const lp_rational_interval_t* I) {
  if (I->is_point) {
    if (!N->is_point) {
      rational_destruct(&N->b);
    }
    rational_neg(&N->a, &I->a);
    N->b_open = N->a_open = 0;
    N->is_point = 1;
    return;
  }

  if (N->is_point) {
    rational_construct(&N->b);
    N->is_point = 0;
  }

  // -[a, b] = [-b, -a]
  // (doing swap in case I and N are the same)
  rational_neg(&N->a, &I->a);
  rational_neg(&N->b, &I->b);
  N->a_open = I->a_open;
  N->b_open = I->b_open;

  rational_swap(&N->a, &N->b);
  size_t tmp = N->a_open;
  N->a_open = N->b_open;
  N->b_open = tmp;
}

void dyadic_interval_neg(lp_dyadic_interval_t* N, const lp_dyadic_interval_t* I) {
  if (I->is_point) {
    if (!N->is_point) {
      dyadic_rational_destruct(&N->b);
    }
    dyadic_rational_neg(&N->a, &I->a);
    N->b_open = N->a_open = 0;
    N->is_point = 1;
    return;
  }

  if (N->is_point) {
    dyadic_rational_construct(&N->b);
    N->is_point = 0;
  }

  // -[a, b] = [-b, -a]
  // (doing swap in case I and N are the same)
  dyadic_rational_neg(&N->a, &I->a);
  dyadic_rational_neg(&N->b, &I->b);
  N->a_open = I->a_open;
  N->b_open = I->b_open;

  dyadic_rational_swap(&N->a, &N->b);
  size_t tmp = N->a_open;
  N->a_open = N->b_open;
  N->b_open = tmp;
}

void rational_interval_sub(lp_rational_interval_t* S, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {
  lp_rational_interval_t neg;
  lp_rational_interval_construct_copy(&neg, I2);
  rational_interval_neg(&neg, &neg);
  rational_interval_add(S, I1, &neg);
  lp_rational_interval_destruct(&neg);
}

void dyadic_interval_sub(lp_dyadic_interval_t* S, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  lp_dyadic_interval_t neg;
  lp_dyadic_interval_construct_copy(&neg, I2);
  dyadic_interval_neg(&neg, &neg);
  dyadic_interval_add(S, I1, &neg);
  lp_dyadic_interval_destruct(&neg);
}

void rational_interval_mul(lp_rational_interval_t* P, const lp_rational_interval_t* I1, const lp_rational_interval_t* I2) {
  if (I1->is_point) {
    if (I2->is_point) {
      // Just multiply the points
      rational_mul(&P->a, &I1->a, &I2->a);
      if (!P->is_point) {
        rational_destruct(&P->b);
        P->is_point = 1;
      }
      P->a_open = P->b_open = 0;
    } else {
      // I1 is point
      // I2 is not point

      // Depending on the sign of a, we might have to flip
      int a_sgn = rational_sgn(&I1->a);
      if (a_sgn == 0) {
        // It's just 0
        if (!P->is_point) {
          rational_destruct(&P->b);
          P->is_point = 1;
        }
        P->a_open = P->b_open = 0;
        rational_assign_int(&P->a, 0, 1);
      } else if (a_sgn > 0) {
        lp_rational_interval_t result;
        rational_construct(&result.a);
        rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        rational_mul(&result.a, &I1->a, &I2->a);
        rational_mul(&result.b, &I1->a, &I2->b);
        lp_rational_interval_swap(&result, P);
        lp_rational_interval_destruct(&result);
      } else {
        lp_rational_interval_t result;
        rational_construct(&result.a);
        rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        rational_mul(&result.b, &I1->a, &I2->a);
        rational_mul(&result.a, &I1->a, &I2->b);
        lp_rational_interval_swap(&result, P);
        lp_rational_interval_destruct(&result);
      }
    }
  } else if (I2->is_point) {
    rational_interval_mul(P, I2, I1);
  } else {
    if (P->is_point) {
      rational_construct(&P->b);
      P->is_point = 0;
    }

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //

    lp_rational_interval_t result;
    lp_rational_interval_construct_zero(&result);

    lp_rational_t tmp;
    rational_construct(&tmp);

    // I1.a x I2.a
    rational_mul(&result.a, &I1->a, &I2->a);
    rational_construct_copy(&result.b, &result.a);
    result.a_open = result.b_open = I1->a_open || I2->a_open;
    result.is_point = 0;

    // I1.a x I2.b
    int tmp_open = I1->a_open || I2->b_open;
    rational_mul(&tmp, &I1->a, &I2->b);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    tmp_open = I1->b_open || I2->a_open;
    rational_mul(&tmp, &I1->b, &I2->a);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    tmp_open = I1->b_open || I2->b_open;
    rational_mul(&tmp, &I1->b, &I2->b);
    if (rational_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (rational_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    lp_rational_interval_swap(&result, P);
    lp_rational_interval_destruct(&result);
    rational_destruct(&tmp);
  }
}

void dyadic_interval_mul(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  if (I1->is_point) {
    if (I2->is_point) {
      // Just multiply the points
      dyadic_rational_mul(&P->a, &I1->a, &I2->a);
      if (!P->is_point) {
        dyadic_rational_destruct(&P->b);
        P->is_point = 1;
      }
      P->a_open = P->b_open = 0;
    } else {
      // Depending on the sign of a, we might have to flip
      int a_sgn = dyadic_rational_sgn(&I1->a);
      if (a_sgn == 0) {
        // It's just 0
        if (!P->is_point) {
          dyadic_rational_destruct(&P->b);
          P->is_point = 1;
        }
        P->a_open = P->b_open = 0;
        dyadic_rational_assign_int(&P->a, 0, 1);
      } else if (a_sgn > 0) {
        // Regular multiplication
        lp_dyadic_interval_t result;
        dyadic_rational_construct(&result.a);
        dyadic_rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->a_open;
        result.b_open = I2->b_open;
        dyadic_rational_mul(&result.a, &I1->a, &I2->a);
        dyadic_rational_mul(&result.b, &I1->a, &I2->b);
        lp_dyadic_interval_swap(&result, P);
        lp_dyadic_interval_destruct(&result);
      } else {
        lp_dyadic_interval_t result;
        dyadic_rational_construct(&result.a);
        dyadic_rational_construct(&result.b);
        result.is_point = 0;
        result.a_open = I2->b_open;
        result.b_open = I2->a_open;
        dyadic_rational_mul(&result.a, &I1->a, &I2->b);
        dyadic_rational_mul(&result.b, &I1->a, &I2->a);
        lp_dyadic_interval_swap(&result, P);
        lp_dyadic_interval_destruct(&result);
      }
    }
  } else if (I2->is_point) {
    dyadic_interval_mul(P, I2, I1);
  } else {
    if (P->is_point) {
      dyadic_rational_construct(&P->b);
      P->is_point = 0;
    }

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //

    lp_dyadic_interval_t result;
    lp_dyadic_interval_construct_zero(&result);

    lp_dyadic_rational_t tmp;
    dyadic_rational_construct(&tmp);

    // I1.a x I2.a
    dyadic_rational_mul(&result.a, &I1->a, &I2->a);
    dyadic_rational_construct_copy(&result.b, &result.a); // was not constructed (is 0)
    result.a_open = result.b_open = I1->a_open || I2->a_open;
    result.is_point = 0;

    // I1.a x I2.b
    int tmp_open = I1->a_open || I2->b_open;
    dyadic_rational_mul(&tmp, &I1->a, &I2->b);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    tmp_open = I1->b_open || I2->a_open;
    dyadic_rational_mul(&tmp, &I1->b, &I2->a);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    tmp_open = I1->b_open || I2->b_open;
    dyadic_rational_mul(&tmp, &I1->b, &I2->b);
    if (dyadic_interval_endpoint_lt(&tmp, tmp_open, &result.a, result.a_open)) {
      dyadic_rational_swap(&tmp, &result.a);
      result.a_open = tmp_open;
    } else if (dyadic_interval_endpoint_lt(&result.b, result.b_open, &tmp, tmp_open)) {
      dyadic_rational_swap(&tmp, &result.b);
      result.b_open = tmp_open;
    }

    lp_dyadic_interval_swap(&result, P);
    lp_dyadic_interval_destruct(&result);
    dyadic_rational_destruct(&tmp);
  }
}

static
int lp_value_mul_approx(const lp_value_t* v1, const lp_value_t* v2, lp_value_t* lb, lp_value_t* ub) {

  assert(v1 != lb && v1 != ub);
  assert(v2 != lb && v2 != ub);
  assert(lb != ub);

  int is_point = 1;
  lp_integer_t mul_int;
  lp_rational_t mul_rat;
  lp_dyadic_rational_t mul_dy;
  lp_dyadic_interval_t mul_dy_interval;

  if (v1->type == v2->type) {
    switch (v1->type) {
    case LP_VALUE_MINUS_INFINITY:
    case LP_VALUE_PLUS_INFINITY:
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_PLUS_INFINITY, 0);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_PLUS_INFINITY, 0);
      }
      break;
    case LP_VALUE_INTEGER:
      lp_integer_construct(&mul_int);
      lp_integer_mul(lp_Z, &mul_int, &v1->value.z, &v2->value.z);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_INTEGER, &mul_int);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_INTEGER, &mul_int);
      }
      lp_integer_destruct(&mul_int);
      break;
    case LP_VALUE_RATIONAL:
      lp_rational_construct(&mul_rat);
      lp_rational_mul(&mul_rat, &v1->value.q, &v2->value.q);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_RATIONAL, &mul_rat);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_RATIONAL, &mul_rat);
      }
      lp_rational_destruct(&mul_rat);
      break;
    case LP_VALUE_DYADIC_RATIONAL:
      lp_dyadic_rational_construct(&mul_dy);
      lp_dyadic_rational_mul(&mul_dy, &v1->value.dy_q, &v2->value.dy_q);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &mul_dy);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &mul_dy);
      }
      lp_dyadic_rational_destruct(&mul_dy);
      break;
    case LP_VALUE_ALGEBRAIC:
      if (v1->value.a.I.is_point && v2->value.a.I.is_point) {
        lp_dyadic_rational_construct(&mul_dy);
        lp_dyadic_rational_mul(&mul_dy, &v1->value.a.I.a, &v2->value.a.I.a);
        if (lb) {
          lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &mul_dy);
        }
        if (ub) {
          lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &mul_dy);
        }
        lp_dyadic_rational_destruct(&mul_dy);
      } else {
        lp_dyadic_interval_construct_zero(&mul_dy_interval);
        dyadic_interval_mul(&mul_dy_interval, &v1->value.a.I, &v2->value.a.I);
        if (lb) {
          lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &mul_dy_interval.a);
        }
        if (ub) {
          lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &mul_dy_interval.b);
        }
        is_point = 0;
        lp_dyadic_interval_destruct(&mul_dy_interval);
      }
      break;
    case LP_VALUE_NONE:
      assert(0);
      break;
    }
  } else {

    lp_value_t v1_tmp, v2_tmp;
    const lp_value_t* v1_to_use;
    const lp_value_t* v2_to_use;

    // Cast to the same type
    int cast = lp_value_to_same_type(v1, v2, &v1_tmp, &v2_tmp, &v1_to_use, &v2_to_use);
    if (cast) {
      // Cast successful, just add as usual
      lp_value_mul_approx(v1_to_use, v2_to_use, lb, ub);
      // If a fesh value was used, delete it
      if (v1_to_use != v1) {
        lp_value_destruct(&v1_tmp);
      }
      if (v2_to_use != v2) {
        lp_value_destruct(&v2_tmp);
      }
    } else {
      // Cast failed, check the cases
      int v1_sgn = lp_value_sgn(v1);
      int v2_sgn = lp_value_sgn(v2);

      if (v1_sgn == 0 || v2_sgn == 0) {
        if (lb) {
          lp_value_assign_zero(lb);
        }
        if (ub) {
          lp_value_assign_zero(ub);
        }
        return 1;
      } else {
        if (lp_value_is_infinity(v1) || lp_value_is_infinity(v2)) {
          int sgn = v1_sgn*v2_sgn;
          lp_value_type_t result = LP_VALUE_NONE;
          if (sgn > 0) {
            result = LP_VALUE_PLUS_INFINITY;
          } else if (sgn < 0) {
            result = LP_VALUE_MINUS_INFINITY;
          } else {
            assert(0);
          }
          if (lb) {
            lp_value_assign_raw(lb, result, 0);
          }
          if (ub) {
            lp_value_assign_raw(ub, result, 0);
          }
          return 1;
        }
      }

       // OK, now deal with algebraic numbers, assume v1 is algebraic
      if (v2->type == LP_VALUE_ALGEBRAIC) {
        const lp_value_t* tmp = v1; v1 = v2; v2 = tmp;
      }

      assert(v1->type == LP_VALUE_ALGEBRAIC);
      assert(v2->type != LP_VALUE_MINUS_INFINITY);
      assert(v2->type != LP_VALUE_PLUS_INFINITY);
      assert(v2->type != LP_VALUE_ALGEBRAIC);

      // v1 is of the form (dy1, dy2), v2 is a proper non-algebraic value
      // we add them into (dy1 * v2, dy2 * v2). if v2 is negative, we swap the bounds

      lp_value_t v1_lb, v1_ub;
      const lp_dyadic_rational_t* a = &v1->value.a.I.a;
      const lp_dyadic_rational_t* b = v1->value.a.I.is_point ? a : &v1->value.a.I.b;
      lp_value_construct(&v1_lb, LP_VALUE_DYADIC_RATIONAL, a);
      lp_value_construct(&v1_ub, LP_VALUE_DYADIC_RATIONAL, b);
      lp_value_mul_approx(&v1_lb, v2, lb, 0);
      lp_value_mul_approx(&v1_ub, v2, ub, 0);
      if (v2_sgn < 0) {
        lp_value_swap(lb, ub);
      }
      is_point = 0;
      lp_value_destruct(&v1_lb);
      lp_value_destruct(&v1_ub);
    }
  }

  return is_point;
}

void lp_interval_mul(lp_interval_t* mul, const lp_interval_t* I1, const lp_interval_t* I2) {

  lp_interval_t result;
  lp_interval_construct_full(&result);

  if (I1->is_point) {
    if (I2->is_point) {
      // Just multiply the points
      result.is_point = lp_value_mul_approx(&I1->a, &I2->a, &result.a, &result.b);
      if (result.is_point) {
        lp_value_destruct(&result.b);
      }
      result.a_open = result.b_open = !result.is_point;
    } else {
      // Depending on the sign of a, we might have to flip
      int a_sgn = lp_value_sgn(&I1->a);
      if (a_sgn == 0) {
        // It's just 0
        lp_interval_destruct(&result);
        lp_interval_construct_zero(&result);
      } else if (a_sgn > 0) {
        // Regular multiplication
        int a_is_point = lp_value_mul_approx(&I1->a, &I2->a, &result.a, 0);
        int b_is_point = lp_value_mul_approx(&I1->a, &I2->b, 0, &result.b);
        result.is_point = 0;
        result.a_open = I2->a_open || !a_is_point;
        result.b_open = I2->b_open || !b_is_point;
      } else {
        int a_is_point = lp_value_mul_approx(&I1->a, &I2->b, &result.a, 0);
        int b_is_point = lp_value_mul_approx(&I1->a, &I2->a, 0, &result.b);
        result.is_point = 0;
        result.a_open = I2->b_open || !a_is_point;
        result.b_open = I2->a_open || !b_is_point;
      }
    }
  } else if (I2->is_point) {
    return lp_interval_mul(mul, I2, I1);
  } else {

    assert(!I1->is_point && !I2->is_point);

    //
    // I1 x I2 = { x*y | x in I1, y in I2 }
    //         = { x*y | I1.a < x < I1.b, I2.a < y < I2.b }
    //

    lp_value_t tmp_lb, tmp_ub;
    lp_value_construct_zero(&tmp_lb);
    lp_value_construct_zero(&tmp_ub);

    int mul_is_point = 0;

    // I1.a x I2.a
    mul_is_point = lp_value_mul_approx(&I1->a, &I2->a, &result.a, &result.b);
    result.a_open = result.b_open = I1->a_open || I2->a_open || !mul_is_point;

    // I1.a x I2.b
    mul_is_point = lp_value_mul_approx(&I1->a, &I2->b, &tmp_lb, &tmp_ub);
    int tmp_open = I1->a_open || I2->b_open || !mul_is_point;
    if (lp_interval_endpoint_lt(&tmp_lb, tmp_open, &result.a, result.a_open)) {
      lp_value_swap(&tmp_lb, &result.a);
      result.a_open = tmp_open;
    }
    if (lp_interval_endpoint_lt(&result.b, result.b_open, &tmp_ub, tmp_open)) {
      lp_value_swap(&tmp_ub, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.a
    mul_is_point = lp_value_mul_approx(&I1->b, &I2->a, &tmp_lb, &tmp_ub);
    tmp_open = I1->b_open || I2->a_open || !mul_is_point;
    if (lp_interval_endpoint_lt(&tmp_lb, tmp_open, &result.a, result.a_open)) {
      lp_value_swap(&tmp_lb, &result.a);
      result.a_open = tmp_open;
    }
    if (lp_interval_endpoint_lt(&result.b, result.b_open, &tmp_ub, tmp_open)) {
      lp_value_swap(&tmp_ub, &result.b);
      result.b_open = tmp_open;
    }

    // I1.b x I2.b
    mul_is_point = lp_value_mul_approx(&I1->b, &I2->b, &tmp_lb, &tmp_ub);
    tmp_open = I1->b_open || I2->b_open || !mul_is_point;
    if (lp_interval_endpoint_lt(&tmp_lb, tmp_open, &result.a, result.a_open)) {
      lp_value_swap(&tmp_lb, &result.a);
      result.a_open = tmp_open;
    }
    if (lp_interval_endpoint_lt(&result.b, result.b_open, &tmp_ub, tmp_open)) {
      lp_value_swap(&tmp_ub, &result.b);
      result.b_open = tmp_open;
    }

    // Final check: if an endpoint is 0, it must have come from another zero
    // endpoint. If that endpoint is closed, then the resulting endpoint must
    // be closed too.
    int sgn_a = lp_value_sgn(&result.a);
    if (sgn_a == 0) {
      int c1_a = (lp_value_sgn(&I1->a) == 0) && !I1->a_open;
      int c1_b = (lp_value_sgn(&I1->b) == 0) && !I1->b_open;
      int c2_a = (lp_value_sgn(&I2->a) == 0) && !I2->a_open;
      int c2_b = (lp_value_sgn(&I2->b) == 0) && !I2->b_open;
      if (c1_a || c1_b || c2_a || c2_b) {
        result.a_open = 0;
      }
    }
    int sgn_b = lp_value_sgn(&result.b);
    if (sgn_b == 0) {
      int c1_a = (lp_value_sgn(&I1->a) == 0) && !I1->a_open;
      int c1_b = (lp_value_sgn(&I1->b) == 0) && !I1->b_open;
      int c2_a = (lp_value_sgn(&I2->a) == 0) && !I2->a_open;
      int c2_b = (lp_value_sgn(&I2->b) == 0) && !I2->b_open;
      if (c1_a || c1_b || c2_a || c2_b) {
        result.b_open = 0;
      }
    }

    lp_value_destruct(&tmp_lb);
    lp_value_destruct(&tmp_ub);
  }

  // Put into result
  lp_interval_swap(mul, &result);
  lp_interval_destruct(&result);
}

static
int lp_value_pow_approx(const lp_value_t* v, unsigned n, lp_value_t* lb, lp_value_t* ub) {

  assert(n > 0);
  assert(v != lb);
  assert(v != ub);
  assert(ub != lb);

  lp_integer_t pow_int;
  lp_rational_t pow_rat;
  lp_dyadic_rational_t pow_dy;
  lp_dyadic_interval_t pow_dy_interval;

  int is_point = 1;

  switch (v->type) {
  case LP_VALUE_INTEGER:
    lp_integer_construct(&pow_int);
    integer_pow(lp_Z, &pow_int, &v->value.z, n);
    if (lb) {
      lp_value_assign_raw(lb, LP_VALUE_INTEGER, &pow_int);
    }
    if (ub) {
      lp_value_assign_raw(ub, LP_VALUE_INTEGER, &pow_int);
    }
    lp_integer_destruct(&pow_int);
    break;
  case LP_VALUE_RATIONAL:
    lp_rational_construct(&pow_rat);
    rational_pow(&pow_rat, &v->value.q, n);
    if (lb) {
      lp_value_assign_raw(lb, LP_VALUE_RATIONAL, &pow_rat);
    }
    if (ub) {
      lp_value_assign_raw(ub, LP_VALUE_RATIONAL, &pow_rat);
    }
    lp_rational_destruct(&pow_rat);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    lp_dyadic_rational_construct(&pow_dy);
    dyadic_rational_pow(&pow_dy, &v->value.dy_q, n);
    if (lb) {
      lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &pow_dy);
    }
    if (ub) {
      lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &pow_dy);
    }
    lp_dyadic_rational_destruct(&pow_dy);
    break;
  case LP_VALUE_ALGEBRAIC:
    if (v->value.a.I.is_point) {
      lp_dyadic_rational_construct(&pow_dy);
      dyadic_rational_pow(&pow_dy, &v->value.a.I.a, n);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &pow_dy);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &pow_dy);
      }
      lp_dyadic_rational_destruct(&pow_dy);
    } else {
      lp_dyadic_interval_construct_zero(&pow_dy_interval);
      dyadic_interval_pow(&pow_dy_interval, &v->value.a.I, n);
      if (lb) {
        lp_value_assign_raw(lb, LP_VALUE_DYADIC_RATIONAL, &pow_dy_interval.a);
      }
      if (ub) {
        lp_value_assign_raw(ub, LP_VALUE_DYADIC_RATIONAL, &pow_dy_interval.b);
      }
      lp_dyadic_interval_destruct(&pow_dy_interval);
      is_point = 0;
    }
    break;
  case LP_VALUE_MINUS_INFINITY:
  case LP_VALUE_PLUS_INFINITY: {
    int sgn = (n % 2) ? lp_value_sgn(v) : 1;
    lp_value_type_t type = sgn > 0 ? LP_VALUE_PLUS_INFINITY : LP_VALUE_MINUS_INFINITY;
    if (lb) {
      lp_value_assign_raw(lb, type, 0);
    }
    if (ub) {
      lp_value_assign_raw(ub, type, 0);
    }
    break;
  }
  case LP_VALUE_NONE:
    assert(0);
  }

  return is_point;
}

int lp_interval_sgn(const lp_interval_t* I) {
  int a_sgn = lp_value_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = lp_value_sgn(&I->b);

  if (a_sgn < 0 && b_sgn > 0) {
    // Definitively contains 0
    return 0;
  }
  if (a_sgn == 0) {
    if (!I->a_open) {
      // Left contains zero
      return 0;
    } else {
      return 1;
    }
  }
  if (b_sgn == 0) {
    if (!I->b_open) {
      // Right contains zero
      return 0;
    } else {
      return -1;
    }
  }

  // Doesn't contain zero
  if (a_sgn < 0) {
    return -1;
  }

  assert(b_sgn > 0);
  return 1;
}

void lp_interval_pow(lp_interval_t* pow, const lp_interval_t* I, unsigned n) {

  lp_interval_t result;
  lp_interval_construct_full(&result);

  if (n == 0) {
    // I^0 = [1]
    lp_value_t one;
    lp_value_construct_int(&one, 1);
    lp_interval_destruct(&result);
    lp_interval_construct_point(&result, &one);
    lp_value_destruct(&one);
  } else if (I->is_point) {
    // Plain power
    result.is_point = lp_value_pow_approx(&I->a, n, &result.a, &result.b);
    if (result.is_point) {
      lp_value_destruct(&result.b);
    }
    result.a_open = result.b_open = !result.is_point;
  } else {
    if (n % 2) {
      // For odd powers we are monotonic, i.e. [a, b]^n = [a^n, b^n]
      int a_point = lp_value_pow_approx(&I->a, n, &result.a, 0);
      int b_point = lp_value_pow_approx(&I->b, n, 0, &result.b);
      result.a_open = I->a_open || !a_point;
      result.b_open = I->b_open || !b_point;
    } else {
      // Even powers depend on whether 0 is in the interval
      int sgn = lp_interval_sgn(I);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        int a_point = lp_value_pow_approx(&I->a, n, 0, &result.a);
        int b_point = lp_value_pow_approx(&I->b, n, 0, &result.b);
        if (lp_interval_endpoint_lt(&result.b, I->b_open, &result.a, I->a_open)) {
          lp_value_swap(&result.b, &result.a);
          result.b_open = I->a_open || !a_point;
        } else {
          result.b_open = I->b_open || !b_point;
        }
        lp_value_assign_zero(&result.a);
        result.a_open = 0;
      } else if (sgn > 0) {
        // P = I^n
        int a_point = lp_value_pow_approx(&I->a, n, &result.a, 0);
        int b_point = lp_value_pow_approx(&I->b, n, 0, &result.b);
        result.a_open = I->a_open || !a_point;
        result.b_open = I->b_open || !b_point;
      } else {
        // P = I^n, but swappeed
        int a_point = lp_value_pow_approx(&I->a, n, 0, &result.b);
        int b_point = lp_value_pow_approx(&I->b, n, &result.a, 0);
        result.a_open = I->b_open || !b_point;
        result.b_open = I->a_open || !a_point;
      }
    }
  }

  lp_interval_swap(pow, &result);
  lp_interval_destruct(&result);
}

void rational_interval_pow(lp_rational_interval_t* P, const lp_rational_interval_t* I, unsigned n) {
  if (n == 0) {
    // I^0 = [1]
    if (!P->is_point) {
      P->is_point = 1;
      rational_destruct(&P->b);
    }
    rational_assign_int(&P->a, 1, 1);
    P->a_open = 0;
    P->b_open = 0;
  } else if (I->is_point) {
    // Plain power
    if (!P->is_point) {
      rational_destruct(&P->b);
      P->is_point = 1;
      P->a_open = P->b_open = 0;
    }
    rational_pow(&P->a, &I->a, n);
  } else {
    if (P->is_point) {
      P->is_point = 0;
      rational_construct(&P->b);
    }
    if (n % 2) {
      // For odd powers we are monotonic, i.e. [a, b]^n = [a^n, b^n]
      P->a_open = I->a_open;
      P->b_open = I->b_open;
      rational_pow(&P->a, &I->a, n);
      rational_pow(&P->b, &I->b, n);
    } else {
      // Even powers depend on whether 0 is in the interval
      int sgn = lp_rational_interval_sgn(I);
      rational_pow(&P->a, &I->a, n);
      rational_pow(&P->b, &I->b, n);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        if (rational_interval_endpoint_lt(&P->b, I->b_open, &P->a, I->a_open)) {
          rational_swap(&P->b, &P->a);
          P->b_open = I->a_open;
        } else {
          P->b_open = I->b_open;
        }
        rational_assign_int(&P->a, 0, 1);
        P->a_open = 0;
      } else if (sgn > 0) {
        // P = I^n
        P->a_open = I->a_open;
        P->b_open = I->b_open;
      } else {
        // negative turns positive, so we flip
        rational_swap(&P->a, &P->b);
        P->a_open = I->b_open;
        P->b_open = I->a_open;
      }
    }
  }
}

void dyadic_interval_pow(lp_dyadic_interval_t* P, const lp_dyadic_interval_t* I, unsigned n) {
  if (n == 0) {
    // I^0 = [1]
    if (!P->is_point) {
      P->is_point = 1;
      dyadic_rational_destruct(&P->b);
    }
    dyadic_rational_assign_int(&P->a, 1, 1);
    P->a_open = 0;
    P->b_open = 0;
  } else if (I->is_point) {
    // Plain power
    if (!P->is_point) {
      dyadic_rational_destruct(&P->b);
      P->is_point = 1;
      P->a_open = P->b_open = 0;
    }
    dyadic_rational_pow(&P->a, &I->a, n);
  } else {
    if (P->is_point) {
      P->is_point = 0;
      dyadic_rational_construct(&P->b);
    }
    if (n % 2) {
      // For odd powers we are monotonic, i.e. [a, b]^n = [a^n, b^n]
      P->a_open = I->a_open;
      P->b_open = I->b_open;
      dyadic_rational_pow(&P->a, &I->a, n);
      dyadic_rational_pow(&P->b, &I->b, n);
    } else {
      // Even powers depend on whether 0 is in the interval
      int sgn = lp_dyadic_interval_sgn(I);
      dyadic_rational_pow(&P->a, &I->a, n);
      dyadic_rational_pow(&P->b, &I->b, n);
      if (sgn == 0) {
        // P = [0, max(a, b)^n]
        if (dyadic_interval_endpoint_lt(&P->b, I->b_open, &P->a, I->a_open)) {
          dyadic_rational_swap(&P->b, &P->a);
          P->b_open = I->a_open;
        } else {
          P->b_open = I->b_open;
        }
        dyadic_rational_assign_int(&P->a, 0, 1);
        P->a_open = 0;
      } else if (sgn > 0) {
        // P = I^n
        P->a_open = I->a_open;
        P->b_open = I->b_open;
      } else {
        // negative turns positive, so we flip
        dyadic_rational_swap(&P->a, &P->b);
        P->a_open = I->b_open;
        P->b_open = I->a_open;
      }
    }
  }
}
