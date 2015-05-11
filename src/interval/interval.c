/*
 * interval.c
 *
 *  Created on: Jan 25, 2014
 *      Author: dejan
 */

#include <interval.h>

#include "number/integer.h"
#include "number/rational.h"
#include "number/dyadic_rational.h"

#include <assert.h>
#include <limits.h>

void lp_interval_construct(lp_interval_t* I,
    const lp_rational_t* a, int a_open,
    const lp_rational_t* b, int b_open)
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

void lp_dyadic_interval_construct(lp_dyadic_interval_t* I,
    const lp_dyadic_rational_t* a, int a_open,
    const lp_dyadic_rational_t* b, int b_open)
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

void lp_interval_construct_point(lp_interval_t* I, const lp_rational_t* a)
{
  rational_construct_copy(&I->a, a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_construct_point(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
  dyadic_rational_construct_copy(&I->a, q);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_interval_construct_zero(lp_interval_t* I) {
  rational_construct(&I->a);
  I->a_open = I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_construct_zero(lp_dyadic_interval_t* I) {
  dyadic_rational_construct(&I->a);
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_interval_construct_copy(lp_interval_t* I, const lp_interval_t* from) {
  rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_dyadic_interval_construct_copy(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from) {
  dyadic_rational_construct_copy(&I->a, &from->a);
  if (!from->is_point) {
    dyadic_rational_construct_copy(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_interval_construct_from_dyadic(lp_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open) {
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

void lp_interval_construct_from_dyadic_interval(lp_interval_t* I, const lp_dyadic_interval_t* from) {
  rational_construct_from_dyadic(&I->a, &from->a);
  if (!from->is_point) {
    rational_construct_from_dyadic(&I->b, &from->b);
  }
  I->a_open = from->a_open;
  I->b_open = from->b_open;
  I->is_point = from->is_point;
}

void lp_dyadic_interval_construct_from_dyadic(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open, const lp_dyadic_rational_t* b, int b_open) {
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

void lp_interval_construct_from_int(lp_interval_t* I,
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

void lp_dyadic_interval_construct_from_int(lp_dyadic_interval_t* I,
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
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_interval_construct_from_integer(lp_interval_t* I,
    const lp_integer_t* a, int a_open,
    const lp_integer_t* b, int b_open)
{
  int cmp = integer_cmp(lp_Z, a, b);
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

void lp_dyadic_interval_construct_from_integer(lp_dyadic_interval_t* I,
    const lp_integer_t* a, int a_open,
    const lp_integer_t* b, int b_open)
{
  int cmp = integer_cmp(lp_Z, a, b);
  assert(cmp <= 0);
  dyadic_rational_construct_from_integer(&I->a, a);
  if (cmp != 0) {
    dyadic_rational_construct_from_integer(&I->b, b);
    I->is_point = 0;
    I->a_open = a_open;
    I->b_open = b_open;
  } else {
    assert(!a_open && !b_open);
    I->is_point = 1;
    I->a_open = I->b_open = 0;
  }
}

void lp_interval_destruct(lp_interval_t* I) {
  rational_destruct(&I->a);
  if (!I->is_point) {
    rational_destruct(&I->b);
  }
}

void lp_dyadic_interval_destruct(lp_dyadic_interval_t* I) {
  dyadic_rational_destruct(&I->a);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
}

void lp_interval_assign(lp_interval_t* I, const lp_interval_t* from) {
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

void lp_dyadic_interval_assign(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* from) {
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

void lp_interval_swap(lp_interval_t* I1, lp_interval_t* I2) {
  lp_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

void lp_dyadic_interval_swap(lp_dyadic_interval_t* I1, lp_dyadic_interval_t* I2) {
  lp_dyadic_interval_t tmp = *I1;
  *I1 = *I2;
  *I2 = tmp;
}

int lp_interval_sgn(const lp_interval_t* I) {
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

int lp_dyadic_interval_sgn(const lp_dyadic_interval_t* I) {
  int a_sgn = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return a_sgn;
  }
  int b_sgn = dyadic_rational_sgn(&I->b);

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

int lp_dyadic_interval_cmp_integer(const lp_dyadic_interval_t* I, const lp_rational_t* z) {

  if (I->is_point) {
    return dyadic_rational_cmp_integer(&I->a, z);
  }

  // I = [a, b]

  int cmp_lower = dyadic_rational_cmp_integer(&I->a, z);
  if (cmp_lower > 0) {
    // a > z => [a, b] > z
    return 1;
  }
  if (cmp_lower == 0) {
    if (I->a_open) {
      // a == z => (a, b] > z
      return 1;
    } else {
      // a == z => [a, b] == z
      return 0;
    }
  }

  int cmp_upper = dyadic_rational_cmp_integer(&I->b, z);
  if (cmp_upper < 0) {
    // [a, b] < z
    return -1;
  }
  if (cmp_upper == 0) {
    if (I->b_open) {
      // [a, b) < z
      return -1;
    } else {
      // [a, b] == z
      return 0;
    }
  }

  // It's inside, return 0
  return 0;
}

int lp_dyadic_interval_cmp_dyadic_rational(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {

}

int lp_dyadic_interval_cmp_rational(const lp_dyadic_interval_t* I, const lp_rational_t* q) {

}


int lp_interval_contains(const lp_interval_t* I, const lp_rational_t* q) {
  int cmp_a = rational_cmp(&I->a, q);
  if (I->is_point) {
    return cmp_a == 0;
  }
  if (I->a_open && !(cmp_a < 0)) return 0;
  if (!I->a_open && !(cmp_a <= 0)) return 0;
  int cmp_b = rational_cmp(q, &I->b);
  if (I->b_open && !(cmp_b < 0)) return 0;
  if (!I->b_open && !(cmp_b <= 0)) return 0;
  return 1;
}

int lp_dyadic_interval_contains(const lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
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

int lp_interval_contains_zero(const lp_interval_t* I) {
  int cmp_a = rational_sgn(&I->a);
  if (I->is_point) {
    return cmp_a == 0;
  }
  if (I->a_open && !(cmp_a < 0)) return 0;
  if (!I->a_open && !(cmp_a <= 0)) return 0;
  int cmp_b = rational_sgn(&I->b);
  if (I->b_open && !(cmp_b < 0)) return 0;
  if (!I->b_open && !(cmp_b <= 0)) return 0;
  return 1;
}

int lp_dyadic_interval_contains_zero(const lp_dyadic_interval_t* I) {
  int cmp_a = dyadic_rational_sgn(&I->a);
  if (I->is_point) {
    return cmp_a == 0;
  }
  if (I->a_open && !(cmp_a < 0)) return 0;
  if (!I->a_open && !(cmp_a <= 0)) return 0;
  int cmp_b = dyadic_rational_sgn(&I->b);
  if (I->b_open && !(cmp_b < 0)) return 0;
  if (!I->b_open && !(cmp_b <= 0)) return 0;
  return 1;
}


void lp_dyadic_interval_construct_from_split(lp_dyadic_interval_t* I_left, lp_dyadic_interval_t* I_right, const lp_dyadic_interval_t* I, int left_open, int right_open) {
  assert(!I->is_point);
  lp_dyadic_rational_t m;
  dyadic_rational_construct(&m);
  dyadic_rational_add(&m, &I->a, &I->b);
  dyadic_rational_div_2exp(&m, &m, 1);
  lp_dyadic_interval_construct(I_left, &I->a, I->a_open, &m, left_open);
  lp_dyadic_interval_construct(I_right, &m, right_open, &I->b, I->b_open);
  dyadic_rational_destruct(&m);
}

void lp_dyadic_interval_construct_intersection(lp_dyadic_interval_t* I, const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  if (I1->is_point) {
    assert(lp_dyadic_interval_contains(I2, &I1->a));
    lp_dyadic_interval_construct_copy(I, I1);
  } else if (I2->is_point) {
    assert(lp_dyadic_interval_contains(I1, &I2->a));
    lp_dyadic_interval_construct_copy(I, I2);
  } else {
    // (   [  )   ]
    int cmp_a = dyadic_rational_cmp(&I1->a, &I2->a);
    const lp_dyadic_rational_t* max_a = cmp_a < 0 ? &I2->a : &I1->a;
    int a_open;
    if (cmp_a == 0) {
      a_open = I1->a_open || I2->a_open;
    } else if (cmp_a < 0) {
      a_open = I2->a_open;
    } else {
      a_open = I1->a_open;
    }

    int cmp_b = dyadic_rational_cmp(&I1->b, &I2->b);
    const lp_dyadic_rational_t* min_b = cmp_b < 0 ? &I1->b : &I2->b;
    int b_open;
    if (cmp_b == 0) {
      b_open = I1->b_open || I2->b_open;
    } else if (cmp_b < 0) {
      b_open = I1->b_open;
    } else {
      b_open = I2->b_open;
    }

    lp_dyadic_interval_construct(I, max_a, a_open, min_b, b_open);
  }
}

void lp_dyadic_interval_collapse_to(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* q) {
  dyadic_rational_assign(&I->a, q);
  if (!I->is_point) {
    dyadic_rational_destruct(&I->b);
  }
  I->a_open = 0;
  I->b_open = 0;
  I->is_point = 1;
}

void lp_dyadic_interval_set_a(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* a, int a_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(a, &I->b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->a, a);
    I->a_open = a_open;
  } else {
    assert(!a_open && !I->b_open);
    lp_dyadic_interval_collapse_to(I, a);
  }
}

void lp_dyadic_interval_set_b(lp_dyadic_interval_t* I, const lp_dyadic_rational_t* b, int b_open) {
  assert(!I->is_point);
  int cmp = dyadic_rational_cmp(&I->a, b);
  assert(cmp <= 0);
  if (cmp != 0) {
    dyadic_rational_assign(&I->b, b);
    I->b_open = b_open;
  } else {
    assert(!I->a_open && !b_open);
    lp_dyadic_interval_collapse_to(I, b);
  }
}

int lp_dyadic_interval_equals(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
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

int lp_dyadic_interval_disjunct(const lp_dyadic_interval_t* I1, const lp_dyadic_interval_t* I2) {
  if (I1->is_point) {
    return !lp_dyadic_interval_contains(I2, &I1->a);
  }
  if (I2->is_point) {
    return !lp_dyadic_interval_contains(I1, &I2->a);
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

void lp_dyadic_interval_scale(lp_dyadic_interval_t* I, int n) {
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
int lp_interval_print(const lp_interval_t* I, FILE* out) {
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

int lp_dyadic_interval_print(const lp_dyadic_interval_t* I, FILE* out) {
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
