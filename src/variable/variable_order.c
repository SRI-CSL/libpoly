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

#include <variable_db.h>
#include <variable_list.h>
#include <variable_order.h>

#include "variable/variable_order.h"

#include <assert.h>
#include <stdlib.h>

/**
 * A simple variable order that orders variable based on a given list, and
 * order the rest of the variables based on their variable id.
 */
struct lp_variable_order_struct {
  /** Reference count */
  size_t ref_count;
  /** The actual order */
  lp_variable_list_t list;
  /** Special top variable */
  lp_variable_t top;
  /** Special bottom variable */
  lp_variable_t bot;
};

void lp_variable_order_construct(lp_variable_order_t* var_order) {
  // No-one pointing yet
  var_order->ref_count = 0;
  // The list
  lp_variable_list_construct(&var_order->list);
  var_order->bot = lp_variable_null;
  var_order->top = lp_variable_null;
}

void lp_variable_order_reverse(lp_variable_order_t* var_order) {
  size_t size = var_order->list.list_size;
  if (size > 1) {
    size_t first, last;
    lp_variable_t* vars = var_order->list.list;
    for (first = 0, last = size-1; first < last; ++ first, --last) {
      lp_variable_t tmp = vars[first];
      vars[first] = vars[last];
      vars[last] = tmp;
    }
  }
}

void lp_variable_order_destruct(lp_variable_order_t* var_order) {
  lp_variable_list_destruct(&var_order->list);
}

void lp_variable_order_attach(lp_variable_order_t* var_order) {
  lp_variable_order_t* self = (lp_variable_order_t*) var_order;
  self->ref_count ++;
}

void lp_variable_order_detach(lp_variable_order_t* var_order) {
  lp_variable_order_t* self = (lp_variable_order_t*) var_order;
  assert(self->ref_count > 0);
  self->ref_count --;
  if (self->ref_count == 0) {
    lp_variable_order_destruct(self);
    free(var_order);
  }
}

lp_variable_order_t* lp_variable_order_new(void) {
  lp_variable_order_t* var_order = malloc(sizeof(lp_variable_order_t));
  lp_variable_order_construct(var_order);
  lp_variable_order_attach((lp_variable_order_t*) var_order);
  return (lp_variable_order_t*) var_order;
}

int lp_variable_order_cmp(const lp_variable_order_t* var_order, lp_variable_t x, lp_variable_t y) {
  const lp_variable_order_t* self = (lp_variable_order_t*) var_order;

  if (x == y) {
    return 0;
  }

  // if bot is smaller than anything
  if (x == var_order->bot) { return -1; }
  if (y == var_order->bot) { return 1;  }
  if (x == var_order->top) { return 1;  }
  if (y == var_order->top) { return -1; }

  // Compare indices
  int x_index = lp_variable_list_index(&self->list, x);
  int y_index = lp_variable_list_index(&self->list, y);

  int cmp = 0;
  if (x_index == y_index) {
    // Indices same, just compare the variables
    cmp = ((int) x) - ((int) y);
  } else {
    // If a variable doesn't have an index, it's bigger
    if (x_index == -1) {
      cmp = 1;
    } else if (y_index == -1) {
      cmp = -1;
    } else {
      cmp = x_index - y_index;
    }
  }

  return cmp;
}

size_t lp_variable_order_size(const lp_variable_order_t* var_order) {
  return lp_variable_list_size(&var_order->list);
}

void lp_variable_order_clear(lp_variable_order_t* var_order) {
  while (lp_variable_list_size(&var_order->list)) {
    lp_variable_list_pop(&var_order->list);
  }
}

void lp_variable_order_push(lp_variable_order_t* var_order, lp_variable_t var) {
  lp_variable_list_push(&var_order->list, var);
}

void lp_variable_order_pop(lp_variable_order_t* var_order) {
  lp_variable_list_pop(&var_order->list);
}

lp_variable_t lp_variable_order_top(const lp_variable_order_t* var_order) {
  return lp_variable_list_top(&var_order->list);
}

int lp_variable_order_print(const lp_variable_order_t* var_order, const lp_variable_db_t* var_db, FILE* out) {
  size_t i;
  int ret = 0;
  ret += fprintf(out, "[");
  for (i = 0; i < var_order->list.list_size; ++ i) {
    if (i) {
      ret += fprintf(out, ", ");
    }
    ret += fprintf(out, "%s", lp_variable_db_get_name(var_db, var_order->list.list[i]));
  }
  ret += fprintf(out, "]");
  return ret;
}

#if _XOPEN_SOURCE >= 700 || _POSIX_C_SOURCE >= 200809L
char* lp_variable_order_to_string(const lp_variable_order_t* var_order, const lp_variable_db_t* var_db) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_variable_order_print(var_order, var_db, f);
  fclose(f);
  return str;
}
#endif

int lp_variable_order_contains(lp_variable_order_t* var_order, lp_variable_t x) {
  return lp_variable_list_index(&var_order->list, x) != -1;
}

void lp_variable_order_make_top(lp_variable_order_t* var_order, lp_variable_t var) {
  assert(var_order->top == lp_variable_null || var == lp_variable_null);
  var_order->top = var;
}

void lp_variable_order_make_bot(lp_variable_order_t* var_order, lp_variable_t var) {
  assert(var_order->bot == lp_variable_null || var == lp_variable_null);
  var_order->bot = var;
}

const lp_variable_list_t* lp_variable_order_get_list(const lp_variable_order_t* var_order) {
  return &var_order->list;
}
