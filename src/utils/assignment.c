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

#include <assignment.h>
#include <variable_db.h>
#include <value.h>
#include <interval.h>

#include "polynomial/polynomial.h"
#include "number/value.h"

#include <stdlib.h>

#define DEFAULT_ASSIGNMENT_SIZE 100

static
void lp_assignment_ensure_size(lp_assignment_t* m, size_t size) {
  if (size > m->size) {
    m->values = realloc(m->values, sizeof(lp_value_t)*size);
    size_t i;
    for (i = m->size; i < size; ++ i) {
      lp_value_construct(m->values + i, LP_VALUE_NONE, 0);
    }
    m->size = size;
  }
}

void lp_assignment_construct(lp_assignment_t* m, const lp_variable_db_t* var_db) {
  m->size = 0;
  m->values = 0;
  m->var_db = var_db;
  lp_variable_db_attach((lp_variable_db_t*)var_db);
  lp_assignment_ensure_size(m, DEFAULT_ASSIGNMENT_SIZE);
}

lp_assignment_t* lp_assignment_new(const lp_variable_db_t* var_db) {
  lp_assignment_t* new = malloc(sizeof(lp_assignment_t));
  lp_assignment_construct(new, var_db);
  return new;
}

void lp_assignment_destruct(lp_assignment_t* m) {
  if (m->values) {
    size_t i;
    for (i = 0; i < m->size; ++ i) {
      lp_value_destruct(m->values + i);
    }
    free(m->values);
  }
  lp_variable_db_detach((lp_variable_db_t*)m->var_db);
}

void lp_assignment_delete(lp_assignment_t* m) {
  lp_assignment_destruct(m);
  free(m);
}

int lp_assignment_print(const lp_assignment_t* m, FILE* out) {
  size_t i, j, ret = 0;
  ret += fprintf(out, "[");
  for (i = 0, j = 0; i < m->size; ++ i) {
    if (m->values[i].type != LP_VALUE_NONE) {
      if (j ++) {
        ret += fprintf(out, ", ");
      }
      ret += fprintf(out, "%s -> ", lp_variable_db_get_name(m->var_db, i));
      ret += lp_value_print(m->values + i, out);
    }
  }
  ret += fprintf(out, "]");
  return ret;
}

char* lp_assignment_to_string(const lp_assignment_t* m) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_assignment_print(m, f);
  fclose(f);
  return str;
}

void lp_assignment_set_value(lp_assignment_t* m, lp_variable_t x, const lp_value_t* value) {
  if (value) {
    lp_assignment_ensure_size(m, x + 1);
    lp_value_destruct(m->values + x);
    lp_value_construct_copy(m->values + x, value);
  } else {
    if (m->size > x) {
      if ((m->values + x)->type != LP_VALUE_NONE) {
        lp_value_destruct(m->values + x);
        lp_value_construct(m->values + x, LP_VALUE_NONE, 0);
      }
    }
  }
}

const lp_value_t* lp_assignment_get_value(const lp_assignment_t* m, lp_variable_t x) {
  if (x < m->size) {
    return m->values + x;
  } else {
    return lp_value_none();
  }
}

void lp_assignment_get_value_approx(const lp_assignment_t* m, lp_variable_t x, lp_rational_interval_t* approx) {
  assert(lp_assignment_get_value(m, x)->type != LP_VALUE_NONE);
  const lp_value_t* x_value = lp_assignment_get_value(m, x);
  lp_value_approx(x_value, approx);
}

int lp_assignment_sgn(const lp_assignment_t* m, const lp_polynomial_t* A) {
  return lp_polynomial_sgn(A, m);
}

static
void lp_interval_assignment_ensure_size(lp_interval_assignment_t* m_const, size_t size) {
  lp_interval_assignment_t* m = (lp_interval_assignment_t*) m_const;
  if (size > m->size) {
    m->intervals = realloc(m->intervals, sizeof(lp_interval_t)*size);
    m->timestamps = realloc(m->timestamps, sizeof(size_t)*size);
    size_t i;
    for (i = m->size; i < size; ++ i) {
      lp_interval_construct_full(m->intervals + i);
      m->timestamps[i] = m->timestamp;
    }
    m->size = size;
  }
}

void lp_interval_assignment_construct(lp_interval_assignment_t* m, const lp_variable_db_t* var_db) {
  m->size = 0;
  m->intervals = 0;
  m->var_db = var_db;
  m->timestamp = 1;
  m->timestamps = 0;
  lp_variable_db_attach((lp_variable_db_t*)var_db);
  lp_interval_assignment_ensure_size(m, DEFAULT_ASSIGNMENT_SIZE);
}

lp_interval_assignment_t* lp_interval_assignment_new(const lp_variable_db_t* var_db) {
  lp_interval_assignment_t* new = malloc(sizeof(lp_interval_assignment_t));
  lp_interval_assignment_construct(new, var_db);
  return new;
}

void lp_interval_assignment_destruct(lp_interval_assignment_t* m) {
  if (m->intervals) {
    size_t i;
    for (i = 0; i < m->size; ++ i) {
      lp_interval_destruct(m->intervals + i);
    }
    free(m->intervals);
    free(m->timestamps);
  }
  lp_variable_db_detach((lp_variable_db_t*)m->var_db);
}

void lp_interval_assignment_delete(lp_interval_assignment_t* m) {
  lp_interval_assignment_destruct(m);
  free(m);
}

int lp_interval_assignment_print(const lp_interval_assignment_t* m, FILE* out) {
  size_t i, j, ret = 0;
  ret += fprintf(out, "[");
  for (i = 0, j = 0; i < m->size; ++ i) {
    if (m->timestamps[i] < m->timestamp) {
      continue;
    }
    if (j ++) {
      ret += fprintf(out, ", ");
    }
    ret += fprintf(out, "%s -> ", lp_variable_db_get_name(m->var_db, i));
    ret += lp_interval_print(m->intervals + i, out);
  }
  ret += fprintf(out, "]");
  return ret;
}

char* lp_interval_assignment_to_string(const lp_interval_assignment_t* m) {
  char* str = 0;
  size_t size = 0;
  FILE* f = open_memstream(&str, &size);
  lp_interval_assignment_print(m, f);
  fclose(f);
  return str;
}

void lp_interval_assignment_set_interval(lp_interval_assignment_t* m, lp_variable_t x, const lp_interval_t* value) {
  lp_interval_assignment_ensure_size(m, x + 1);
  lp_interval_destruct(m->intervals + x);
  if (value) {
    lp_interval_construct_copy(m->intervals + x, value);
  } else {
    lp_interval_construct_full(m->intervals + x);
  }
  m->timestamps[x] = m->timestamp;
}

const lp_interval_t* lp_interval_assignment_get_interval(const lp_interval_assignment_t* m, lp_variable_t x) {
  if (x >= m->size) {
    return lp_interval_full();
  }
  if (m->timestamps[x] == m->timestamp) {
    return m->intervals + x;
  } else {
    return lp_interval_full();
  }
}

void lp_interval_assignment_reset(lp_interval_assignment_t* m) {
  m->timestamp ++;
}
