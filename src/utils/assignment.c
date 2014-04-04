/*
 * assignment_internal.c
 *
 *  Created on: Mar 12, 2014
 *      Author: dejan
 */

#include "utils/assignment.h"

#include "number/value.h"
#include "polynomial/polynomial.h"

#include <malloc.h>

#define DEFAULT_ASSIGNMENT_SIZE 100
#define DEFAULT_ASSIGNMENT_SIZE 100

static
void assignment_ensure_size(const assignment_t* m_const, size_t size) {
  assignment_t* m = (assignment_t*) m_const;
  if (size >= m->size) {
    m->values = realloc(m->values, sizeof(value_t)*size);
    size_t i;
    for (i = m->size; i < size; ++ i) {
      value_construct(m->values + i, VALUE_NONE, 0);
    }
    m->size = size;
  }
}

void assignment_construct(assignment_t* m, const variable_db_t* var_db) {
  m->size = 0;
  m->values = 0;
  m->var_db = var_db;
  variable_db_ops.attach((variable_db_t*)var_db);
  assignment_ensure_size(m, DEFAULT_ASSIGNMENT_SIZE);
}

assignment_t* assignment_new(const variable_db_t* var_db) {
  assignment_t* new = malloc(sizeof(assignment_t));
  assignment_construct(new, var_db);
  return new;
}

void assignment_destruct(assignment_t* m) {
  if (m->values) {
    size_t i;
    for (i = 0; i < m->size; ++ i) {
      value_ops.destruct(m->values + i);
    }
  }
  variable_db_ops.detach((variable_db_t*)m->var_db);
}

int assignment_print(const assignment_t* m, FILE* out) {
  int i, j, ret = 0;
  ret += fprintf(out, "[");
  for (i = 0, j = 0; i < ret; ++ i) {
    if (m->values[i].type != VALUE_NONE) {
      if (j ++) {
        ret += fprintf(out, ", ");
      }
      ret += fprintf(out, "%s -> ", variable_db_ops.get_name(m->var_db, i));
      ret += value_ops.print(m->values + i, out);
    }
  }
  ret += fprintf(out, "]");
  return ret;
}

char* assignment_to_string(const assignment_t* m) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  assignment_print(m, f);
  fclose(f);
  return str;
}


void assignment_set_value(assignment_t* m, variable_t x, const value_t* value) {
  if (value) {
    assignment_ensure_size(m, x + 1);
    value_destruct(m->values + x);
    value_construct_copy(m->values + x, value);
  } else {
    if (m->size > x) {
      if (m->values->type != VALUE_NONE) {
        value_destruct(m->values + x);
        value_construct(m->values + x, VALUE_NONE, 0);
      }
    }
  }
}

const value_t* assignment_get_value(const assignment_t* m, variable_t x) {
  assignment_ensure_size(m, x + 1);
  return m->values + x;
}

void assignment_get_value_approx(const assignment_t* m, variable_t x, interval_t* approx) {
  assert(assignment_get_value(m, x)->type != VALUE_NONE);
  const value_t* x_value = assignment_get_value(m, x);
  value_approx(x_value, approx);
}

int assignment_sgn(const assignment_t* m, const polynomial_t* A) {
  return polynomial_sgn(A, m);
}

//
// API Construction
//

const assignment_ops_t assignment_ops = {
    assignment_construct,
    assignment_new,
    assignment_destruct,
    assignment_print,
    assignment_to_string,
    assignment_set_value,
    assignment_get_value,
    assignment_get_value_approx,
    assignment_sgn
};
