/*
 * variable_order.c
 *
 *  Created on: Feb 7, 2014
 *      Author: dejan
 */

#include "variable/variable_order.h"

#include <assert.h>
#include <malloc.h>

void variable_order_simple_construct(variable_order_simple_t* var_order) {
  // No-one pointing yet
  var_order->ref_count = 0;
  // The operations
  var_order->ops = &variable_order_simple_ops;
  // The list
  variable_list_ops.construct(&var_order->list);
}

void variable_order_simple_destruct(variable_order_simple_t* var_order) {
  variable_list_ops.destruct(&var_order->list);
}

void variable_order_simple_attach(variable_order_t* var_order) {
  variable_order_simple_t* self = (variable_order_simple_t*) var_order;
  self->ref_count ++;
}

void variable_order_simple_detach(variable_order_t* var_order) {
  variable_order_simple_t* self = (variable_order_simple_t*) var_order;
  assert(self->ref_count > 0);
  self->ref_count --;
  if (self->ref_count == 0) {
    variable_order_simple_destruct(self);
    free(var_order);
  }
}

variable_order_t* variable_order_simple_new(void) {
  variable_order_simple_t* var_order = malloc(sizeof(variable_order_simple_t));
  variable_order_simple_construct(var_order);
  variable_order_simple_attach((variable_order_t*) var_order);
  return (variable_order_t*) var_order;
}

int variable_order_simple_cmp(const variable_order_t* var_order, variable_t x, variable_t y) {
  const variable_order_simple_t* self = (variable_order_simple_t*) var_order;

  int x_index = variable_list_ops.index(&self->list, x);
  int y_index = variable_list_ops.index(&self->list, y);

  if (x_index == y_index) {
    return ((int) x) - ((int) y);
  } else {
    return x_index - y_index;
  }
}

size_t variable_order_simple_size(const variable_order_simple_t* var_order) {
  return variable_list_ops.size(&var_order->list);
}

void variable_order_simple_clear(variable_order_simple_t* var_order) {
  while (variable_list_ops.size(&var_order->list)) {
    variable_list_ops.pop(&var_order->list);
  }
}

void variable_order_simple_push(variable_order_simple_t* var_order, variable_t var) {
  variable_list_ops.push(&var_order->list, var);
}

void variable_order_simple_pop(variable_order_simple_t* var_order) {
  variable_list_ops.pop(&var_order->list);
}

int variable_order_simple_print(const variable_order_simple_t* var_order, const variable_db_t* var_db, FILE* out) {
  size_t i;
  int ret = 0;
  ret += fprintf(out, "[");
  for (i = 0; i < var_order->list.list_size; ++ i) {
    if (i) {
      ret += fprintf(out, ", ");
    }
    ret += fprintf(out, "%s", variable_db_ops.get_name(var_db, var_order->list.list[i]));
  }
  ret += fprintf(out, "]");
  return ret;
}

char* variable_order_simple_to_string(const variable_order_simple_t* var_order, const variable_db_t* var_db) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  variable_order_simple_print(var_order, var_db, f);
  fclose(f);
  return str;
}

int variable_order_simple_contains(variable_order_simple_t* var_order, variable_t x) {
  return variable_list_ops.index(&var_order->list, x) != -1;
}

variable_order_simple_ops_t variable_order_simple_ops = {
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
