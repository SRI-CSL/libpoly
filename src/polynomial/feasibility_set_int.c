/**
 * Copyright 2024, SRI International.
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

#include <feasibility_set_int.h>
#include "polynomial/feasibility_set_int.h"

#include <stdlib.h>
#include <assert.h>

#include "utils/u_memstream.h"

static
int int_cmp(const void *i1, const void *i2) {
  return lp_integer_cmp(lp_Z, i1, i2);
}

static
void lp_feasibility_set_int_construct(lp_int_ring_t *K, lp_feasibility_set_int_t* s, size_t size, bool inverted) {
  s->K = K;
  lp_int_ring_attach(s->K);
  s->size = size;
  s->inverted = inverted;
  s->elements = size ? malloc(s->size * sizeof(lp_integer_t)) : NULL;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_new_empty(lp_int_ring_t *K) {
  lp_feasibility_set_int_t *s = malloc(sizeof(lp_feasibility_set_int_t));
  lp_feasibility_set_int_construct(K, s, 0, false);
  return s;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_new_full(lp_int_ring_t *K) {
  lp_feasibility_set_int_t *s = malloc(sizeof(lp_feasibility_set_int_t));
  lp_feasibility_set_int_construct(K, s, 0, true);
  return s;
}

static
void lp_feasibility_set_int_construct_from_integer(lp_feasibility_set_int_t *set, lp_int_ring_t *K, const lp_integer_t *elements, size_t size, bool inverted) {
  lp_feasibility_set_int_construct(K, set, size, inverted);
  for (size_t i = 0; i < size; ++ i) {
    lp_integer_construct_copy(K, set->elements + i, elements + i);
  }
  set->size = size;
  qsort(set->elements, size, sizeof(lp_integer_t), int_cmp);
}

lp_feasibility_set_int_t* lp_feasibility_set_int_new_from_integer(lp_int_ring_t *K, const lp_integer_t* elements, size_t elements_size, bool inverted) {
  lp_feasibility_set_int_t *result = malloc(sizeof(lp_feasibility_set_int_t));
  lp_feasibility_set_int_construct_from_integer(result, K, elements, elements_size, inverted);
  return result;
}

static
void lp_feasibility_set_int_construct_copy(lp_feasibility_set_int_t *set, const lp_feasibility_set_int_t* from) {
  // we assume that from is valid, thus no sorting
  lp_feasibility_set_int_construct(from->K, set, from->size, from->inverted);
  for (size_t i = 0; i < from->size; ++ i) {
    lp_integer_construct_copy(lp_Z, set->elements + i, from->elements + i);
  }
  set->size = from->size;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_new_copy(const lp_feasibility_set_int_t* set) {
  lp_feasibility_set_int_t *result = malloc(sizeof(lp_feasibility_set_int_t));
  lp_feasibility_set_int_construct_copy(result, set);
  return result;
}

static
void lp_feasibility_set_int_destruct(lp_feasibility_set_int_t* set) {
  lp_int_ring_detach(set->K);
  free(set->elements);
}

void lp_feasibility_set_int_delete(lp_feasibility_set_int_t* set) {
  assert(set);
  lp_feasibility_set_int_destruct(set);
  free(set);
}

void lp_feasibility_set_int_assign(lp_feasibility_set_int_t* set, const lp_feasibility_set_int_t* from) {
  if (set != from) {
    lp_feasibility_set_int_destruct(set);
    lp_feasibility_set_int_construct_copy(set, from);
  }
}

void lp_feasibility_set_int_swap(lp_feasibility_set_int_t* s1, lp_feasibility_set_int_t* s2) {
  lp_feasibility_set_int_t tmp = *s1;
  *s1 = *s2;
  *s2 = tmp;
}

int lp_feasibility_set_int_is_empty(const lp_feasibility_set_int_t* set) {
  if (!set->inverted && set->size == 0) return 1;
  if (set->inverted && lp_integer_cmp_int(lp_Z, &set->K->M, set->size) == 0) return 1;
  return 0;
}

int lp_feasibility_set_int_is_full(const lp_feasibility_set_int_t* set) {
  if (set->inverted && set->size == 0) return 1;
  if (!set->inverted && lp_integer_cmp_int(lp_Z, &set->K->M, set->size) == 0) return 1;
  return 0;
}

void lp_feasibility_set_int_size(const lp_feasibility_set_int_t *set, lp_integer_t *out) {
  lp_integer_assign_int(lp_Z, out, set->size);
  if (set->inverted) {
    lp_integer_sub(lp_Z, out, &set->K->M, out);
  }
}

static
size_t lp_feasibility_set_int_size_approx(const lp_feasibility_set_int_t *set) {
  if (!set->inverted) {
    return set->size;
  } else {
    if (mpz_fits_ulong_p(&set->K->M)) {
      return lp_integer_to_int(&set->K->M) - set->size;
    } else {
      // size can't be big enough for an actual size of 0 or 1
      return ULONG_MAX;
    }
  }
}

int lp_feasibility_set_int_is_point(const lp_feasibility_set_int_t* set) {
  assert(lp_integer_cmp_int(lp_Z, &set->K->M, 2) > 0);
  if (set->inverted && set->size == 1) return 1;
  if (!set->inverted && lp_integer_cmp_int(lp_Z, &set->K->M, set->size + 1) == 0) return 1;
  return 0;
}

/** returns true if value is in elements. Assumes that the elements are sorted. */
static
bool lp_feasibility_set_int_find(const lp_feasibility_set_int_t* set, const lp_integer_t* value, size_t *pos) {
  if (set->size == 0){
    return false;
  }
  assert(set->elements);
  // including l and r
  size_t l = 0, r = set->size - 1;
  while (l <= r) {
    size_t p = (r + l) >> 1;
    int cmp = lp_integer_cmp(lp_Z, set->elements + p, value);
    if (cmp < 0) {
      // continue left
      r = p - 1;
    } else if (cmp > 0) {
      // continue right
      l = p + 1;
    } else {
      // found
      if (pos)
        *pos = p;
      return true;
    }
  }
  return false;
}

int lp_feasibility_set_int_contains(const lp_feasibility_set_int_t* set, const lp_integer_t* value) {
  bool found = lp_feasibility_set_int_find(set, value, NULL);
  return found != set->inverted;
}

void lp_feasibility_set_int_pick_value(const lp_feasibility_set_int_t* set, lp_integer_t* value) {
  assert(!lp_feasibility_set_int_is_empty(set));
  if (!set->inverted) {
    size_t pos = random() % set->size;
    lp_integer_assign(lp_Z, value, set->elements + pos);
  } else {
    assert(false);
    // TODO implement
    // get random element between (incl.) 0 and (|K| - size)
    // find and add number of <= element
  }
}

/** Calculates i1 \cup i2 */
static
lp_feasibility_set_int_status_t ordered_integer_set_union(
        lp_integer_t **result, size_t *result_size,
        const lp_integer_t *i1, size_t i1_size,
        const lp_integer_t *i2, size_t i2_size) {

  size_t max_size = i1_size + i2_size;
  *result_size = 0;
  if (i1_size == 0 && i2_size == 0) {
    return LP_FEASIBILITY_SET_INT_EMPTY;
  }

  *result = malloc(max_size * sizeof(lp_integer_t));

  bool just_i1 = true, just_i2 = true;
  size_t p1 = 0, p2 = 0;
  while (p1 < i1_size && p2 < i2_size) {
    int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
    if (cmp < 0) {
      lp_integer_construct_copy(lp_Z, *result + *result_size, i1 + p1);
      just_i2 = false;
      p1 ++;
    } else if (cmp > 0) {
      lp_integer_construct_copy(lp_Z, *result + *result_size, i2 + p2);
      just_i1 = false;
      p2 ++;
    } else {
      lp_integer_construct_copy(lp_Z, *result + *result_size, i1 + p1);
      p1 ++; p2 ++;
    }
    (*result_size) ++;
  }
  while (p1 < i1_size) {
    lp_integer_construct_copy(lp_Z, *result + *result_size, i1 + p1++);
    (*result_size) ++;
    just_i2 = false;
  }
  while (p2 < i2_size) {
    lp_integer_construct_copy(lp_Z, *result + *result_size, i2 + p2++);
    (*result_size) ++;
    just_i1 = false;
  }
  *result = realloc(*result, *result_size * sizeof(lp_integer_t));

  if (just_i1) {
    return LP_FEASIBILITY_SET_INT_S1;
  } else if (just_i2) {
    return LP_FEASIBILITY_SET_INT_S2;
  } else {
    return LP_FEASIBILITY_SET_INT_NEW;
  }
}

/** Calculates i1 \cap i2 */
static
lp_feasibility_set_int_status_t ordered_integer_set_intersect(
        lp_integer_t **result, size_t *result_size,
        const lp_integer_t *i1, size_t i1_size,
        const lp_integer_t *i2, size_t i2_size) {

  size_t max_size = i1_size < i2_size ? i1_size : i2_size;
  *result_size = 0;
  if (i1_size == 0 || i2_size == 0) {
    return LP_FEASIBILITY_SET_INT_EMPTY;
  }

  *result = malloc(max_size * sizeof(lp_integer_t));

  size_t p1 = 0, p2 = 0;
  bool all_i1 = true, all_i2 = true;
  while (p1 < i1_size && p2 < i2_size) {
    int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
    if (cmp < 0) {
      p2 ++;
      all_i2 = false;
    } else if (cmp > 0) {
      p1 ++;
      all_i1 = false;
    } else {
      lp_integer_construct_copy(lp_Z, *result + *result_size, i1 + p1);
      p1 ++; p2 ++;
      (*result_size) ++;
    }
  }

  *result = realloc(*result, *result_size * sizeof(lp_integer_t));

  if (*result_size == 0) {
    return LP_FEASIBILITY_SET_INT_EMPTY;
  } else if (all_i1) {
    return LP_FEASIBILITY_SET_INT_S1;
  } else if (all_i2) {
    return LP_FEASIBILITY_SET_INT_S2;
  } else {
    return LP_FEASIBILITY_SET_INT_NEW;
  }
}

/** Calculates i1 \setminus i2 */
static
lp_feasibility_set_int_status_t ordered_integer_set_minus(
        lp_integer_t **result, size_t *result_size,
        const lp_integer_t *i1, size_t i1_size,
        const lp_integer_t *i2, size_t i2_size) {

  size_t max_size = i1_size;
  *result_size = 0;
  *result = malloc(max_size * sizeof(lp_integer_t));

  size_t p1 = 0, p2 = 0;
  while (p1 < i1_size) {
    bool found = false;
    while(p2 < i2_size) {
      // as long i2's value is smaller than current i1, we continue
      int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
      if (cmp == 0) {
        found = true;
      }
      if (cmp <= 0) {
        break;
      }
      p2 ++;
    }
    if (!found) {
      lp_integer_construct_copy(lp_Z, *result + *result_size, i1 + p1);
      (*result_size) ++;
    }
    p1 ++;
  }
  *result = realloc(*result, *result_size * sizeof(lp_integer_t));

  if (*result_size == 0) {
    return LP_FEASIBILITY_SET_INT_EMPTY;
  } else if (*result_size == i1_size) {
    return LP_FEASIBILITY_SET_INT_S1;
  } else {
    return LP_FEASIBILITY_SET_INT_NEW;
  }
}

static
void lp_feasibility_set_int_invert(lp_feasibility_set_int_t *set) {
  assert(set->K != lp_Z);
  // don't invert big fields
  assert(lp_integer_cmp_int(lp_Z, &set->K->M, 10000) < 0);

  size_t cnt = lp_integer_to_int(&set->K->M);
  assert(set->size <= cnt);
  cnt -= set->size;

  lp_integer_t *new = malloc(cnt * sizeof(lp_integer_t));
  lp_integer_t *old = set->elements;
  size_t pos_new = 0;
  size_t pos_old = 0;

  long lb = lp_integer_to_int(&set->K->lb);
  long ub = lp_integer_to_int(&set->K->ub);
  for (long val = lb; val <= ub; ++val) {
    assert(pos_new < cnt);
    assert(pos_old >= set->size || lp_integer_cmp_int(lp_Z, old + pos_old, val) >= 0);
    if (pos_old < set->size && lp_integer_cmp_int(lp_Z, old + pos_old, val) == 0) {
      ++ pos_old;
    } else {
      lp_integer_construct_from_int(lp_Z, new + pos_new, val);
      ++ pos_new;
    }
  }

  free(old);
  set->elements = new;
  set->size = cnt;
  set->inverted = !set->inverted;
}

static inline
lp_feasibility_set_int_status_t invert_i1_i2(lp_feasibility_set_int_status_t status) {
  if (status == LP_FEASIBILITY_SET_INT_S1) { return LP_FEASIBILITY_SET_INT_S2; }
  if (status == LP_FEASIBILITY_SET_INT_S2) { return LP_FEASIBILITY_SET_INT_S1; }
  return status;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_intersect(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2) {
  lp_feasibility_set_int_status_t status;
  return lp_feasibility_set_int_intersect_with_status(s1, s2, &status);
}

lp_feasibility_set_int_t* lp_feasibility_set_int_intersect_with_status(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2, lp_feasibility_set_int_status_t * status) {
  assert(lp_int_ring_equal(s1->K, s2->K));

  lp_feasibility_set_int_t *result;

  if (s1->inverted && s2->inverted) {
    result = lp_feasibility_set_int_new_empty(s1->K);
    *status = ordered_integer_set_union(&result->elements, &result->size,
                                  s1->elements, s1->size, s2->elements,
                                  s2->size);
    result->inverted = true;
  } else if (!s1->inverted && !s2->inverted) {
    result = lp_feasibility_set_int_new_empty(s1->K);
    *status = ordered_integer_set_intersect(&result->elements, &result->size,
                                      s1->elements, s1->size, s2->elements,
                                      s2->size);
    result->inverted = false;
  } else if (s1->inverted && !s2->inverted) {
    result = lp_feasibility_set_int_intersect_with_status(s2, s1, status);
    // TODO this gives I2 precedence in case I1 == I2. Maybe introduce a BOTH option.
    *status = invert_i1_i2(*status);
  } else {
    assert(!s1->inverted && s2->inverted);
    if (s1->size > lp_feasibility_set_int_size_approx(s2)) {
      // s2 could be the smaller set, LP_FEASIBILITY_SET_INT_S2 is possible
      // TODO this effort could be saved in case we don't care about the status
      lp_feasibility_set_int_t *tmp = lp_feasibility_set_int_new_copy(s2);
      lp_feasibility_set_int_invert(tmp);
      result = lp_feasibility_set_int_intersect_with_status(s1, tmp, status);
      lp_feasibility_set_int_delete(tmp);
    } else {
      // s1 is the smaller set, LP_FEASIBILITY_SET_INT_S2 is not possible
      result = lp_feasibility_set_int_new_empty(s1->K);
      *status = ordered_integer_set_minus(&result->elements, &result->size,
                                          s1->elements, s1->size,
                                          s2->elements, s2->size);
      result->inverted = false;
    }
  }

  return result;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_union(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2) {
  lp_feasibility_set_int_status_t status;
  return lp_feasibility_set_int_union_with_status(s1, s2, &status);
}

lp_feasibility_set_int_t* lp_feasibility_set_int_union_with_status(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2, lp_feasibility_set_int_status_t* status) {
  assert(lp_int_ring_equal(s1->K, s2->K));

  lp_feasibility_set_int_t *result = lp_feasibility_set_int_new_empty(s1->K);

  if (s1->inverted && s2->inverted) {
    *status = ordered_integer_set_intersect(&result->elements, &result->size,
                                            s1->elements, s1->size,
                                            s2->elements, s2->size);
    result->inverted = true;
  } else if (!s1->inverted && !s2->inverted) {
    *status = ordered_integer_set_union(&result->elements, &result->size,
                                        s1->elements, s1->size,
                                        s2->elements, s2->size);
    result->inverted = false;
  } else if (!s1->inverted && s2->inverted) {
    result = lp_feasibility_set_int_union_with_status(s2, s1, status);
    // TODO this gives I2 precedence in case I1 == I2. Maybe introduce a BOTH option.
    *status = invert_i1_i2(*status);
  } else {
    assert (s1->inverted && !s2->inverted);
    if (s2->size > lp_feasibility_set_int_size_approx(s1)) {
      // s2 is be the bigger set, LP_FEASIBILITY_SET_INT_S2 is possible
      // TODO this effort could be saved in case we don't care about the status
      lp_feasibility_set_int_t *tmp = lp_feasibility_set_int_new_copy(s1);
      lp_feasibility_set_int_invert(tmp);
      result = lp_feasibility_set_int_intersect_with_status(tmp, s2, status);
      lp_feasibility_set_int_delete(tmp);
    } else {
      // s1 is the bigger set, LP_FEASIBILITY_SET_INT_S2 is not possible
      *status = ordered_integer_set_minus(&result->elements, &result->size,
                                          s1->elements, s1->size,
                                          s2->elements, s2->size);
      result->inverted = true;
    }
  }

  return result;
}

void lp_feasibility_set_int_add(lp_feasibility_set_int_t* s, const lp_feasibility_set_int_t* from) {
  lp_feasibility_set_int_t *tmp = lp_feasibility_set_int_union(s, from);
  lp_feasibility_set_int_swap(tmp, s);
  lp_feasibility_set_int_delete(tmp);
}

int lp_feasibility_set_int_print(const lp_feasibility_set_int_t* set, FILE* out) {
  int ret = 0;
  if (set->inverted) {
    ret += fprintf(out, "F \\");
  }
  ret += fprintf(out, "{ ");
  for(size_t i = 0; i < set->size; ++ i) {
    if (i) {
      ret += fprintf(out, ", ");
    }
    ret += lp_integer_print(&set->elements[i], out);
  }
  ret += fprintf(out, " } in ");
  ret += lp_int_ring_print(set->K, out);
  return ret;
}

char* lp_feasibility_set_int_to_string(const lp_feasibility_set_int_t* set) {
  struct u_memstream mem;
  char* str = 0;
  size_t size = 0;
  u_memstream_open(&mem, &str, &size);
  FILE* f = u_memstream_get(&mem);
  lp_feasibility_set_int_print(set, f);
  u_memstream_close(&mem);
  return str;
}
