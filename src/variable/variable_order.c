/*
 * variable_order.c
 *
 *  Created on: Feb 7, 2014
 *      Author: dejan
 */

#include "variable/variable_order.h"

#include <assert.h>
#include <malloc.h>

void variable_order_simple_construct(lp_variable_order_simple_t* var_order) {
  // No-one pointing yet
  var_order->ref_count = 0;
  // The operations
  var_order->ops = &lp_variable_order_simple_ops;
  // The list
  lp_variable_list_construct(&var_order->list);
}

void variable_order_simple_destruct(lp_variable_order_simple_t* var_order) {
  lp_variable_list_destruct(&var_order->list);
}

void variable_order_simple_attach(lp_variable_order_t* var_order) {
  lp_variable_order_simple_t* self = (lp_variable_order_simple_t*) var_order;
  self->ref_count ++;
}

void variable_order_simple_detach(lp_variable_order_t* var_order) {
  lp_variable_order_simple_t* self = (lp_variable_order_simple_t*) var_order;
  assert(self->ref_count > 0);
  self->ref_count --;
  if (self->ref_count == 0) {
    variable_order_simple_destruct(self);
    free(var_order);
  }
}

lp_variable_order_t* variable_order_simple_new(void) {
  lp_variable_order_simple_t* var_order = malloc(sizeof(lp_variable_order_simple_t));
  variable_order_simple_construct(var_order);
  variable_order_simple_attach((lp_variable_order_t*) var_order);
  return (lp_variable_order_t*) var_order;
}

int variable_order_simple_cmp(const lp_variable_order_t* var_order, lp_variable_t x, lp_variable_t y) {
  const lp_variable_order_simple_t* self = (lp_variable_order_simple_t*) var_order;

  int x_index = lp_variable_list_index(&self->list, x);
  int y_index = lp_variable_list_index(&self->list, y);

  if (x_index == y_index) {
    return ((int) x) - ((int) y);
  } else {
    return x_index - y_index;
  }
}

size_t variable_order_simple_size(const lp_variable_order_simple_t* var_order) {
  return lp_variable_list_size(&var_order->list);
}

void variable_order_simple_clear(lp_variable_order_simple_t* var_order) {
  while (lp_variable_list_size(&var_order->list)) {
    lp_variable_list_pop(&var_order->list);
  }
}

void variable_order_simple_push(lp_variable_order_simple_t* var_order, lp_variable_t var) {
  lp_variable_list_push(&var_order->list, var);
}

void variable_order_simple_pop(lp_variable_order_simple_t* var_order) {
  lp_variable_list_pop(&var_order->list);
}

int variable_order_simple_print(const lp_variable_order_simple_t* var_order, const lp_variable_db_t* var_db, FILE* out) {
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

char* variable_order_simple_to_string(const lp_variable_order_simple_t* var_order, const lp_variable_db_t* var_db) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  variable_order_simple_print(var_order, var_db, f);
  fclose(f);
  return str;
}

int variable_order_simple_contains(lp_variable_order_simple_t* var_order, lp_variable_t x) {
  return lp_variable_list_index(&var_order->list, x) != -1;
}

lp_variable_order_simple_ops_t lp_variable_order_simple_ops = {
    {
        variable_order_simple_new,
        variable_order_simple_attach,
        variable_order_simple_detach,
        variable_order_simple_cmp
    },
    variable_order_simple_size,
    variable_order_simple_clear,
    variable_order_simple_contains,
    variable_order_simple_push,
    variable_order_simple_pop,
    variable_order_simple_print,
    variable_order_simple_to_string
};
