/*
 * interval.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include "interval/interval.h"

#include <assert.h>
#include <limits.h>

void interval_construct(interval_t* I,
    const rational_t* a, int a_open,
    const rational_t* b, int b_open)
{
  int cmp = rational_cmp(a, b);
  assert(cmp <= 0);
  rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    rational_construct_copy(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    I->is_point = 1;
    assert(!a_open && !b_open);
    I->a_open = I->b_open = 0;
  }
}

void dyadic_interval_construct(dyadic_interval_t* I,
    const dyadic_rational_t* a, int a_open,
    const dyadic_rational_t* b, int b_open)
{
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_copy(&I->b, b);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void interval_construct_point(interval_t* I, const rational_t* a)
{
  rational_construct_copy(&I->a, a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void dyadic_interval_construct_point(dyadic_interval_t* I, const dyadic_rational_t* q) {
  dyadic_rational_construct_copy(&I->a, q);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void interval_construct_zero(interval_t* I) {
  rational_construct(&I->a);
  I->a_open = I->b_open = 0;
  I->is_point = 1;
}

void dyadic_interval_construct_zero(dyadic_interval_t* I) {
  dyadic_rational_construct(&I->a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void interval_construct_copy(interval_t* I, const interval_t* from) {
  rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void dyadic_interval_construct_copy(dyadic_interval_t* I, const dyadic_interval_t* from) {
  dyadic_rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    dyadic_rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void interval_construct_from_dyadic(interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  rational_construct_from_dyadic(&I->a, a);
  if (cmp != 0) {
    rational_construct_from_dyadic(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void dyadic_interval_construct_from_dyadic(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open, const dyadic_rational_t* b, int b_open) {
  int cmp = dyadic_rational_cmp(a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_copy(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_copy(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void interval_construct_from_int(interval_t* I,
    long a, int a_open,
    long b, int b_open)
{
  assert(a <= b);
  rational_construct_from_int(&I->a, a, 0);
  if (a != b) {
    rational_construct_from_int(&I->b, b, 0);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void dyadic_interval_construct_from_int(dyadic_interval_t* I,
    long a, int a_open,
    long b, int b_open)
{
  assert(a <= b);
  dyadic_rational_construct_from_int(&I->a, a, 0);
  if (a != b) {
    dyadic_rational_construct_from_int(&I->b, b, 0);
    I->a_open = a_open;
    I->b_open = b_open;
    I->is_point = 0;
  } else {
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void interval_construct_from_integer(interval_t* I,
    const integer_t* a, int a_open,
    const integer_t* b, int b_open)
{
  int cmp = integer_cmp(Z, a, b);
  assert(cmp <= 0);
  rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    rational_construct_from_integer(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void dyadic_interval_construct_from_integer(dyadic_interval_t* I,
    const integer_t* a, int a_open,
    const integer_t* b, int b_open)
{
  int cmp = integer_cmp(Z, a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_from_integer(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void interval_destruct(interval_t* I) {
  rational_destruct(&I->a);
  if (!I->is_point) {
    rational_destruct(&I->b);
  }
}

void dyadic_interval_destruct(dyadic_interval_t* I) {
  dyadic_rational_destruct(&I->a);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
}

void interval_assign(interval_t* I, const interval_t* from) {
  if (I != from) {
    if (I->is_point) {
      if (from->is_point) {
        rational_assign(&I->a, &from->a);
      } else {
        rational_assign(&I->a, &from->a);
        rational_construct_copy(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    } else {
      if (from->is_point) {
        rational_assign(&I->a, &from->a);
        rational_destruct(&I->b);
        I->a_open = I->b_open = 0;
        I->is_point = 1;
      } else {
        // both intervals
        rational_assign(&I->a, &from->a);
        rational_assign(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    }
  }
}

void dyadic_interval_assign(dyadic_interval_t* I, const dyadic_interval_t* from) {
  if (I != from) {
    if (I->is_point) {
      if (from->is_point) {
        dyadic_rational_assign(&I->a, &from->a);
      } else {
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_construct_copy(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    } else {
      if (from->is_point) {
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_destruct(&I->b);
        I->a_open = I->b_open = 0;
        I->is_point = 1;
      } else {
        // both intervals
        dyadic_rational_assign(&I->a, &from->a);
        dyadic_rational_assign(&I->b, &from->b);
        I->a_open = from->a_open;
        I->b_open = from->b_open;
        I->is_point = 0;
      }
    }
  }
}

void interval_swap(interval_t* I1, interval_t* I2) {
  interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

void dyadic_interval_swap(dyadic_interval_t* I1, dyadic_interval_t* I2) {
  dyadic_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

int interval_sgn(const interval_t* I) {
  int a_sgn = rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = rational_sgn(&I->b);

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

void dyadic_interval_construct_from_split(dyadic_interval_t* I_left, dyadic_interval_t* I_right, const dyadic_interval_t* I, int left_open, int right_open) {
  assert(!I->is_point);
  dyadic_rational_t m;
  dyadic_rational_construct(&m);
  dyadic_rational_add(&m, &I->a, &I->b);
  dyadic_rational_div_2exp(&m, &m, 1);
  dyadic_interval_construct(I_left, &I->a, I->a_open, &m, left_open);
  dyadic_interval_construct(I_right, &m, right_open, &I->b, I->b_open);
  dyadic_rational_destruct(&m);
}

int dyadic_interval_contains(const dyadic_interval_t* I, const dyadic_rational_t* q) {
  int cmp_a = dyadic_rational_cmp(&I->a, q);
  if (I->is_point) {
    return cmp_a == 0;
  }
  if (I->a_open && !(cmp_a < 0)) return 0;
  if (!I->a_open && !(cmp_a <= 0)) return 0;
  int cmp_b = dyadic_rational_cmp(q, &I->b);
  if (I->b_open && !(cmp_b < 0)) return 0;
  if (!I->b_open && !(cmp_b <= 0)) return 0;
  return 1;
}

void dyadic_interval_construct_intersection(dyadic_interval_t* I, const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point) {
    assert(dyadic_interval_contains(I2, &I1->a));
    dyadic_interval_construct_copy(I, I1);
  } else if (I2->is_point) {
    assert(dyadic_interval_contains(I1, &I2->a));
    dyadic_interval_construct_copy(I, I2);
  } else {
    // (   [  )   ]
    int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
    const dyadic_rational_t* max_a = cmp_a < 0 ? &I2->a : &I1->a;
    int a_open;
    if (cmp_a == 0) {
      a_open = I1->a_open || I2->a_open;
    } else if (cmp_a < 0) {
      a_open = I2->a_open;
    } else {
      a_open = I1->a_open;
    }

    int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
    const dyadic_rational_t* min_b = cmp_b < 0 ? &I1->b : &I2->b;
    int b_open;
    if (cmp_b == 0) {
      b_open = I1->b_open || I2->b_open;
    } else if (cmp_b < 0) {
      b_open = I1->b_open;
    } else {
      b_open = I2->b_open;
    }

    dyadic_interval_construct(I, max_a, a_open, min_b, b_open);
  }
}

int dyadic_interval_sgn(const dyadic_interval_t* I) {
  int a_sgn = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = dyadic_rational_sgn(&I->b);
  if (a_sgn < 0 && b_sgn <= 0) {
    return -1;
  } else if (b_sgn > 0 && a_sgn >= 0) {
    return 1;
  } else {
    // a_sgn <= 0, b_sgn >= 0, contains 0
    return 0;
  }
}

void dyadic_interval_collapse_to(dyadic_interval_t* I, const dyadic_rational_t* q) {
  dyadic_rational_assign(&I->a, q);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void dyadic_interval_set_a(dyadic_interval_t* I, const dyadic_rational_t* a, int a_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(a, &I->b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->a, a);
    I->a_open = a_open;
  } else {
    assert(!a_open && !I->b_open);
    dyadic_interval_collapse_to(I, a);
  }
}

void dyadic_interval_set_b(dyadic_interval_t* I, const dyadic_rational_t* b, int b_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(&I->a, b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->b, b);
    I->b_open = b_open;
  } else {
    assert(!I->a_open && !b_open);
    dyadic_interval_collapse_to(I, b);
  }
}

int dyadic_interval_equals(const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point && !I2->is_point) {
    return 0;
  }
  if (!I1->is_point && I2->is_point) {
    return 0;
  }
  int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
  if (I1->is_point) {
    assert(I2->is_point);
    return cmp_a == 0;
  }
  if (cmp_a != 0 || ((!I1->a_open) != (!I2->a_open))) return 0;
  int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
  if (cmp_b != 0 || ((!I1->b_open) != (!I2->b_open))) return 0;
  return 1;
}

int dyadic_interval_disjunct(const dyadic_interval_t* I1, const dyadic_interval_t* I2) {
  if (I1->is_point) {
    return !dyadic_interval_contains(I2, &I1->a);
  }
  if (I2->is_point) {
    return !dyadic_interval_contains(I1, &I2->a);
  }
  // ( I1 ) ( I2 ) ?
  int cmp1 = dyadic_rational_cmp(&I1->b, &I2->a);
  if (cmp1 < 0) return 1;
  else if (cmp1 == 0 && (I1->b_open || I2->a_open)) return 1;
  // ( I2 ) ( I1 ) ?
  int cmp2 = dyadic_rational_cmp(&I2->b, &I1->a);
  if (cmp2 < 0) return 1;
  else if (cmp2 == 0 && (I2->b_open || I1->a_open)) return 1;
  return 0;
}

void dyadic_interval_scale(dyadic_interval_t* I, int n) {
  assert(!I->is_point);
  if (n > 0) {
    dyadic_rational_mul_2exp(&I->a, &I->a, n);
    dyadic_rational_mul_2exp(&I->b, &I->b, n);
  } else {
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
    dyadic_rational_div_2exp(&I->a, &I->a, -n);
  }
}

/** Prints the interval to the given stream. */
int interval_print(const interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += rational_print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += rational_print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += rational_print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

int dyadic_interval_print(const dyadic_interval_t* I, FILE* out) {
  int ret = 0;
  if (I) {
    if (I->is_point) {
      ret += fprintf(out, "[");
      ret += dyadic_rational_print(&I->a, out);
      ret += fprintf(out, "]");
    } else {
      if (I->a_open) {
        ret += fprintf(out, "(");
      } else {
        ret += fprintf(out, "[");
      }
      ret += dyadic_rational_print(&I->a, out);
      ret += fprintf(out, ", ");
      ret += dyadic_rational_print(&I->b, out);
      if (I->b_open) {
        ret += fprintf(out, ")");
      } else {
        ret += fprintf(out, "]");
      }
    }
  } else {
    ret += fprintf(out, "(-inf, +inf)");
  }
  return ret;
}

const interval_ops_t interval_ops = {
    interval_construct,
    interval_construct_point,
    interval_construct_copy,
    interval_construct_from_dyadic,
    interval_construct_from_int,
    interval_construct_from_integer,
    interval_destruct,
    interval_swap,
    interval_print
};

const dyadic_interval_ops_t dyadic_interval_ops = {
    dyadic_interval_construct,
    dyadic_interval_construct_copy,
    dyadic_interval_construct_from_int,
    dyadic_interval_construct_from_integer,
    dyadic_interval_construct_point,
    dyadic_interval_construct_from_split,
    dyadic_interval_construct_intersection,
    dyadic_interval_destruct,
    dyadic_interval_swap,
    dyadic_interval_collapse_to,
    dyadic_interval_set_a,
    dyadic_interval_set_b,
    dyadic_interval_print,
    dyadic_interval_equals,
    dyadic_interval_contains,
    dyadic_interval_disjunct,
    dyadic_interval_scale
};
