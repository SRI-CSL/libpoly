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

static
int int_cmp(const void *i1, const void *i2) {
  return lp_integer_cmp(lp_Z, i1, i2);
}

static
void lp_feasibility_set_int_ensure_capacity(lp_feasibility_set_int_t* s, size_t capacity) {
  if (capacity && capacity > s->capacity) {
    s->capacity = capacity;
    s->elements = realloc(s->elements, s->capacity * sizeof(lp_integer_t));
  }
}

static
void lp_feasibility_set_int_construct(lp_int_ring_t *K, lp_feasibility_set_int_t* s, size_t size, bool inverted) {
  s->K = K;
  lp_int_ring_attach(s->K);
  s->capacity = 0;
  s->size = size;
  s->inverted = inverted;
  s->elements = NULL;
  lp_feasibility_set_int_ensure_capacity(s, size);
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
size_t ordered_integer_set_union(lp_integer_t **result, const lp_integer_t *i1, size_t i1_size, const lp_integer_t *i2, size_t i2_size) {
  size_t max_size = i1_size + i2_size;
  *result = malloc(max_size * sizeof(lp_integer_t));
  size_t result_size = 0;

  size_t p1 = 0, p2 = 0;
  while (p1 < i1_size && p2 < i2_size) {
    int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
    if (cmp < 0) {
      lp_integer_construct_copy(lp_Z, *result + result_size, i1 + p1);
      p1 ++;
    } else if (cmp > 0) {
      lp_integer_construct_copy(lp_Z, *result + result_size, i2 + p2);
      p2 ++;
    } else {
      lp_integer_construct_copy(lp_Z, *result + result_size, i1 + p1);
      p1 ++; p2 ++;
    }
    result_size ++;
  }
  while (p1 < i1_size) {
    lp_integer_construct_copy(lp_Z, *result + result_size++, i1 + p1++);
  }
  while (p2 < i2_size) {
    lp_integer_construct_copy(lp_Z, *result + result_size++, i2 + p2++);
  }
  *result = realloc(*result, result_size * sizeof(lp_integer_t));
  return result_size;
}

/** Calculates i1 \cap i2 */
static
size_t ordered_integer_set_intersect(lp_integer_t **result, const lp_integer_t *i1, size_t i1_size, const lp_integer_t *i2, size_t i2_size) {
  size_t max_size = i1_size < i2_size ? i1_size : i2_size;
  *result = malloc(max_size * sizeof(lp_integer_t));
  size_t result_size = 0;

  size_t p1 = 0, p2 = 0;
  while (p1 < i1_size && p2 < i2_size) {
    int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
    if (cmp < 0) {
      p2 ++;
    } else if (cmp > 0) {
      p1 ++;
    } else {
      lp_integer_construct_copy(lp_Z, *result + result_size, i1 + p1);
      p1 ++; p2 ++;
      result_size ++;
    }
  }
  *result = realloc(*result, result_size * sizeof(lp_integer_t));
  return result_size;
}

/** Calculates i1 \setminus i2 */
static
size_t ordered_integer_set_minus(lp_integer_t **result, const lp_integer_t *i1, size_t i1_size, const lp_integer_t *i2, size_t i2_size) {
  size_t max_size = i1_size;
  *result = malloc(max_size * sizeof(lp_integer_t));
  size_t result_size = 0;

  size_t p1 = 0, p2 = 0;
  while (p1 < i1_size) {
    bool i2_found = false;
    while(p2 < i2_size) {
      // as long i2's value is smaller than current i1, we continue
      int cmp = lp_integer_cmp(lp_Z, i1 + p1, i2 + p2);
      if (cmp == 0) {
        i2_found = true;
      }
      if (cmp <= 0) {
        break;
      }
      p2 ++;
    }
    if (!i2_found) {
      lp_integer_construct_copy(lp_Z, *result + result_size, i1 + p1);
      result_size ++;
    }
    p1 ++;
  }
  *result = realloc(*result, result_size * sizeof(lp_integer_t));
  return result_size;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_intersect(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2) {
  assert(lp_int_ring_equal(s1->K, s2->K));

  lp_feasibility_set_int_t *result = lp_feasibility_set_int_new_empty(s1->K);
  size_t s;

  if (s1->inverted && s2->inverted) {
    result->inverted = true;
    s = ordered_integer_set_union(&result->elements,
                                  s1->elements,s1->size,
                                  s2->elements, s2->size);
  } else if (s1->inverted && !s2->inverted) {
    result->inverted = false;
    s = ordered_integer_set_minus(&result->elements,
                                  s2->elements, s2->size,
                                  s1->elements, s1->size);
  } else if (!s1->inverted && s2->inverted) {
    result->inverted = false;
    s = ordered_integer_set_minus(&result->elements,
                                  s1->elements, s1->size,
                                  s2->elements, s2->size);
  } else {
    assert(!s1->inverted && !s2->inverted);
    result->inverted = false;
    s = ordered_integer_set_intersect(&result->elements,
                                      s1->elements, s1->size,
                                      s2->elements, s2->size);
  }
  result->capacity = result->size = s;
  return result;
}

lp_feasibility_set_int_t* lp_feasibility_set_int_intersect_with_status(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2, lp_feasibility_set_int_intersect_status_t* status) {

}

void lp_feasibility_set_int_add(lp_feasibility_set_int_t* s, const lp_feasibility_set_int_t* from) {
  lp_feasibility_set_int_t *tmp = lp_feasibility_set_int_union(s, from);
  lp_feasibility_set_int_swap(tmp, s);
  lp_feasibility_set_int_delete(tmp);
}

lp_feasibility_set_int_t* lp_feasibility_set_int_union(const lp_feasibility_set_int_t* s1, const lp_feasibility_set_int_t* s2) {
  assert(lp_int_ring_equal(s1->K, s2->K));

  lp_feasibility_set_int_t *result = lp_feasibility_set_int_new_empty(s1->K);
  size_t s;

  if (s1->inverted && s2->inverted) {
    result->inverted = true;
    s = ordered_integer_set_intersect(&result->elements,
                                      s1->elements,s1->size,
                                      s2->elements, s2->size);
  } else if (s1->inverted && !s2->inverted) {
    result->inverted = true;
    s = ordered_integer_set_minus(&result->elements,
                                  s1->elements, s1->size,
                                  s2->elements, s2->size);
  } else if (!s1->inverted && s2->inverted) {
    result->inverted = true;
    s = ordered_integer_set_minus(&result->elements,
                                  s2->elements, s2->size,
                                  s1->elements, s1->size);
  } else {
    assert(!s1->inverted && !s2->inverted);
    result->inverted = false;
    s = ordered_integer_set_union(&result->elements,
                                  s1->elements, s1->size,
                                  s2->elements, s2->size);
  }
  result->capacity = result->size = s;
  return result;
}
