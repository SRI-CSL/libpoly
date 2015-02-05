/*
 * variable_order.c
 *
 *  Created on: Feb 7, 2014
 *      Author: dejan
 */

#include <variable_order.h>

#include <assert.h>
#include <malloc.h>

/**
 * A simple variable order that orders variable based on a given list, and
 * order the rest of the variables based on their variable id.
 */
struct lp_variable_order_struct {
  /** The operations */
  lp_variable_order_ops_t* ops;
  /** Reference count */
  size_t ref_count;
  /** The actual order */
  lp_variable_list_t list;
};

void lp_variable_order_construct(lp_variable_order_t* var_order) {
  // No-one pointing yet
  var_order->ref_count = 0;
  // The list
  lp_variable_list_construct(&var_order->list);
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

  int x_index = lp_variable_list_index(&self->list, x);
  int y_index = lp_variable_list_index(&self->list, y);

  if (x_index == y_index) {
    return ((int) x) - ((int) y);
  } else {
    return x_index - y_index;
  }
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

char* lp_variable_order_to_string(const lp_variable_order_t* var_order, const lp_variable_db_t* var_db) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_variable_order_print(var_order, var_db, f);
  fclose(f);
  return str;
}

int lp_variable_order_contains(lp_variable_order_t* var_order, lp_variable_t x) {
  return lp_variable_list_index(&var_order->list, x) != -1;
}
