/*
 * lp_polynomial_hash_set.c
 *
 *  Created on: May 27, 2015
 *      Author: dejan
 */

#include <polynomial_hash_set.h>
#include <polynomial.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

/** Default initial size (must be a power of 2) */
#define LP_POLYNOMIAL_HASH_SET_DEFAULT_SIZE 64

/** Resize threshold: the size is doubled when nelems >= size * RESIZE_RATIO */
#define LP_POLYNOMIAL_HASH_SET_RESIZE_RATIO 0.7

void lp_polynomial_hash_set_construct(lp_polynomial_hash_set_t* set) {
  set->data = malloc(LP_POLYNOMIAL_HASH_SET_DEFAULT_SIZE*sizeof(lp_polynomial_t*));
  memset(set->data, 0, LP_POLYNOMIAL_HASH_SET_DEFAULT_SIZE*sizeof(lp_polynomial_t*));
  set->data_size = LP_POLYNOMIAL_HASH_SET_DEFAULT_SIZE;
  set->size = 0;
  set->resize_threshold = LP_POLYNOMIAL_HASH_SET_DEFAULT_SIZE*LP_POLYNOMIAL_HASH_SET_RESIZE_RATIO;
}

void lp_polynomial_hash_set_destruct(lp_polynomial_hash_set_t* set) {
  // Close the set
  lp_polynomial_hash_set_close(set);
  // Remove all the polynomials
  size_t i;
  for (i = 0; i < set->size; ++ i) {
    lp_polynomial_delete(set->data[i]);
  }
  // Free the data
  free(set->data);
  set->data = NULL;
}

static
void lp_polynomial_hash_set_insert_move(lp_polynomial_t** data, size_t mask, lp_polynomial_t* p) {
  size_t i = lp_polynomial_hash(p) & mask;
  while (data[i] != 0) {
    i ++;
    i &= mask;
  }
  data[i] = p;
}

static
int lp_polynomial_hash_set_insert_copy(lp_polynomial_t** data, size_t mask, const lp_polynomial_t* p) {
  size_t i = lp_polynomial_hash(p) & mask;
  while (data[i] != 0) {
    if (lp_polynomial_eq(data[i], p)) { return 0; }
    i ++;
    i &= mask;
  }
  data[i] = lp_polynomial_new_copy(p);
  return 1;
}

static
int lp_polynomial_hash_set_search(lp_polynomial_t** data, size_t mask, const lp_polynomial_t* p) {
  size_t i = lp_polynomial_hash(p) & mask;
  for (;;) {
    if (lp_polynomial_eq(data[i], p)) return 1;
    if (data[i] == 0) return 0;
    i ++;
    i &= mask;
  }
  return 0;
}

/** Double the size of the set. */
static
void lp_polynomial_hash_set_extend(lp_polynomial_hash_set_t* set) {

  size_t old_data_size = set->data_size;
  size_t new_data_size = old_data_size << 1;

  lp_polynomial_t** new_data = malloc(new_data_size*sizeof(lp_polynomial_t*));
  memset(new_data, 0, new_data_size * sizeof(lp_polynomial_t*));

  size_t mask = new_data_size - 1;
  size_t i;
  for (i = 0; i < old_data_size; ++ i) {
    lp_polynomial_t* p = set->data[i];
    if (p != 0) {
      lp_polynomial_hash_set_insert_move(new_data, mask, p);
    }
  }

  free(set->data);
  set->data = new_data;
  set->data_size = new_data_size;
  set->resize_threshold = new_data_size*LP_POLYNOMIAL_HASH_SET_RESIZE_RATIO;
}

int lp_polynomial_hash_set_contains(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p) {
  assert(p);
  return lp_polynomial_hash_set_search(set->data, set->data_size-1, p);
}

int lp_polynomial_hash_set_insert(lp_polynomial_hash_set_t* set, const lp_polynomial_t* p) {
  assert(p);
  assert(set->data_size > set->size);

  int result = lp_polynomial_hash_set_insert_copy(set->data, set->data_size-1, p);
  if (result) {
    set->size ++;
    if (set->size > set->resize_threshold) {
      lp_polynomial_hash_set_extend(set);
    }
  }

  return result;
}

void lp_polynomial_hash_set_close(lp_polynomial_hash_set_t* set) {
  size_t data_size = set->data_size;
  lp_polynomial_t** data = set->data;

  size_t i, j;
  for (i = 0, j = 0; j < data_size; ++ j) {
    lp_polynomial_t* p = data[j];
    if (p != 0) {
      data[i] = p;
      i ++;
    }
  }

  assert(i == set->size && i < data_size);
}

void lp_polynomial_hash_set_clear(lp_polynomial_hash_set_t* set) {
  lp_polynomial_hash_set_destruct(set);
  lp_polynomial_hash_set_construct(set);
}
