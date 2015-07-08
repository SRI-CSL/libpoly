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

#include <variable_list.h>
#include <assert.h>
#include <stdlib.h>
#include <variable_order.h>

#define INITIAL_LIST_SIZE 100
#define INITIAL_MAP_SIZE 100

static
void lp_variable_list_resize(lp_variable_list_t* list, size_t capacity) {
  assert(capacity > list->list_capacity);
  list->list = (lp_variable_t*) realloc(list->list, capacity*sizeof(lp_variable_t));
  list->list_capacity = capacity;
}

static
void lp_variable_map_resize(lp_variable_list_t* list, size_t capacity) {
  assert(capacity > list->var_to_index_map_capacity);
  list->var_to_index_map = (int*) realloc(list->var_to_index_map, capacity*sizeof(int));
  size_t i;
  for (i = list->var_to_index_map_capacity; i < capacity; ++ i) {
    list->var_to_index_map[i] = -1;
  }
  list->var_to_index_map_capacity = capacity;
}

void lp_variable_list_construct(lp_variable_list_t* list) {
  // The list
  list->list = 0;
  list->list_size = 0;
  list->list_capacity = 0;
  lp_variable_list_resize(list, INITIAL_LIST_SIZE);
  // Map to indices
  list->var_to_index_map = 0;
  list->var_to_index_map_capacity = 0;
  lp_variable_map_resize(list, INITIAL_MAP_SIZE);
}

void lp_variable_list_destruct(lp_variable_list_t* list) {
  free(list->list);
  free(list->var_to_index_map);
  list->list = 0;
  list->list_size = 0;
  list->list_capacity = 0;
}

size_t lp_variable_list_size(const lp_variable_list_t* list) {
  return list->list_size;
}

int lp_variable_list_index(const lp_variable_list_t* list, lp_variable_t x) {
  if (x >= list->var_to_index_map_capacity) {
    return -1;
  } else {
    return list->var_to_index_map[x];
  }
}

int lp_variable_list_contains(const lp_variable_list_t* list, lp_variable_t x) {
  return lp_variable_list_index(list, x) != -1;
}

void lp_variable_list_remove(lp_variable_list_t* list, lp_variable_t x) {
  int index = lp_variable_list_index(list, x);
  if (index != -1) {
    list->list[index] = lp_variable_null;
    list->var_to_index_map[x] = -1;
  }
}

void lp_variable_list_copy_into(const lp_variable_list_t* list, lp_variable_t* vars) {
  size_t i;
  for (i = 0; i < list->list_size; ++ i) {
    vars[i] = list->list[i];
  }
}

void lp_variable_list_push(lp_variable_list_t* list, lp_variable_t var) {
  if (list->list_size == list->list_capacity) {
    lp_variable_list_resize(list, list->list_capacity*2);
  }
  if (var >= list->var_to_index_map_capacity) {
    lp_variable_map_resize(list, var+1);
  }
  assert(list->var_to_index_map[var] == -1);
  list->var_to_index_map[var] = list->list_size;
  list->list[list->list_size ++] = var;
}

void lp_variable_list_pop(lp_variable_list_t* list) {
  assert(list->list_size > 0);
  lp_variable_t var = list->list[-- list->list_size];
  list->var_to_index_map[var] = -1;
}

lp_variable_t lp_variable_list_top(const lp_variable_list_t* list) {
  assert(list->list_size > 0);
  return list->list[list->list_size-1];
}

static
const lp_variable_order_t* lp_variable_list_cmp_order = 0;

int lp_variable_list_cmp(const void* x, const void* y) {
  lp_variable_t x_var = *(lp_variable_t*)x;
  lp_variable_t y_var = *(lp_variable_t*)y;
  return lp_variable_order_cmp(lp_variable_list_cmp_order, x_var, y_var);
}

void lp_variable_list_order(lp_variable_list_t* list, const lp_variable_order_t* order) {
  // Compact
  size_t i, to_keep;
  for (i = 0, to_keep = 0; i < list->list_size; ++ i) {
    lp_variable_t x = list->list[i];
    if (x != lp_variable_null) {
      list->list[to_keep ++] = x;
    }
  }
  list->list_size = to_keep;
  // Sort the list
  lp_variable_list_cmp_order = order;
  qsort(list->list, list->list_size,  sizeof(lp_variable_t), lp_variable_list_cmp);
  // Reconstruct indices
  for (i = 0; i < list->list_size; ++ i) {
    list->var_to_index_map[list->list[i]] = i;
  }
}
