/*
 * variable.c
 *
 *Created on: Jan 28, 2014
 *Author: dejan
 */

#include <variable_db.h>

#include "utils/debug_trace.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>

struct lp_variable_db_struct {
  /** Reference count */
  size_t ref_count;
  /** Size of the current variable database */
  size_t size;
  /** The capacity of the database */
  size_t capacity;
  /** Names of individual variables */
  char** variable_names;
};

static
void lp_variable_db_resize(lp_variable_db_t* var_db, size_t capacity) {
  assert(var_db);
  assert(capacity > var_db->capacity);
  var_db->variable_names = realloc(var_db->variable_names, capacity*sizeof(char*));
  var_db->capacity = capacity;
  size_t i;
  for (i = var_db->size; i < capacity; ++ i) {
    var_db->variable_names[i] = 0;
  }
}

#define INITIAL_VAR_DB_SIZE 100

lp_variable_t lp_variable_db_new_variable(lp_variable_db_t* var_db, const char* name) {
  assert(var_db);
  if (var_db->size == var_db->capacity) {
    lp_variable_db_resize(var_db, 2*var_db->capacity);
  }
  var_db->variable_names[var_db->size] = strdup(name);
  return var_db->size ++;
}

void lp_variable_db_add_variable(lp_variable_db_t* var_db, lp_variable_t var, const char* name) {
  assert(var_db);
  while (var >= var_db->capacity) {
    lp_variable_db_resize(var_db, 2*var_db->capacity);
  }
  assert(var_db->variable_names[var] == 0);
  var_db->variable_names[var] = strdup(name);
}

void lp_variable_db_construct(lp_variable_db_t* var_db) {
  assert(var_db);
  var_db->ref_count = 0;
  var_db->size = 0;
  var_db->capacity = 0;
  var_db->variable_names = 0;
  lp_variable_db_resize(var_db, INITIAL_VAR_DB_SIZE);
}

void lp_variable_db_destruct(lp_variable_db_t* var_db) {
  assert(var_db);
  size_t i;
  for (i = 0; i < var_db->size; ++ i) {
    if (var_db->variable_names[i]) {
      free(var_db->variable_names[i]);
    }
  }
}

void lp_variable_db_attach(lp_variable_db_t* var_db) {
  assert(var_db);
  var_db->ref_count ++;
}

void lp_variable_db_detach(lp_variable_db_t* var_db) {
  assert(var_db);
  assert(var_db->ref_count > 0);
  var_db->ref_count --;
  if (var_db->ref_count == 0) {
    lp_variable_db_destruct(var_db);
    free(var_db);
  }
}

lp_variable_db_t* lp_variable_db_new(void) {
  lp_variable_db_t* new = malloc(sizeof(lp_variable_db_t));
  lp_variable_db_construct(new);
  lp_variable_db_attach(new);
  return new;
}

int lp_variable_db_print(const lp_variable_db_t* var_db, FILE* out) {
  assert(var_db);
  int ret = 0;
  size_t i;
  for (i = 0; i < var_db->size; ++ i) {
    if (var_db->variable_names[i]) {
      ret += fprintf(out, "[%zu] = %s\n", i, var_db->variable_names[i]);
    }
  }
  return ret;
}

const char* lp_variable_db_get_name(const lp_variable_db_t* var_db, lp_variable_t var) {
  assert(var_db);
  assert(var < var_db->size);
  return var_db->variable_names[var];
}
