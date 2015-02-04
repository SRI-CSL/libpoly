/*
 * variable.c
 *
 *Created on: Jan 28, 2014
 *Author: dejan
 */

#include <variable.h>

#include "utils/debug_trace.h"

#include <assert.h>
#include <string.h>
#include <malloc.h>

static inline
void variable_db_resize(lp_variable_db_t* var_db, size_t capacity) {
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

lp_variable_t variable_db_new_variable(lp_variable_db_t* var_db, const char* name) {
  assert(var_db);
  if (var_db->size == var_db->capacity) {
    variable_db_resize(var_db, 2*var_db->capacity);
  }
  var_db->variable_names[var_db->size] = strdup(name);
  return var_db->size ++;
}

void variable_db_add_variable(lp_variable_db_t* var_db, lp_variable_t var, const char* name) {
  assert(var_db);
  while (var >= var_db->capacity) {
    variable_db_resize(var_db, 2*var_db->capacity);
  }
  assert(var_db->variable_names[var] == 0);
  var_db->variable_names[var] = strdup(name);
}

void variable_db_construct(lp_variable_db_t* var_db) {
  assert(var_db);
  var_db->ref_count = 0;
  var_db->size = 0;
  var_db->capacity = 0;
  var_db->variable_names = 0;
  variable_db_resize(var_db, INITIAL_VAR_DB_SIZE);
}

void variable_db_destruct(lp_variable_db_t* var_db) {
  assert(var_db);
  size_t i;
  for (i = 0; i < var_db->size; ++ i) {
    if (var_db->variable_names[i]) {
      free(var_db->variable_names[i]);
    }
  }
}

void variable_db_attach(lp_variable_db_t* var_db) {
  assert(var_db);
  var_db->ref_count ++;
}

void variable_db_detach(lp_variable_db_t* var_db) {
  assert(var_db);
  assert(var_db->ref_count > 0);
  var_db->ref_count --;
  if (var_db->ref_count == 0) {
    variable_db_destruct(var_db);
    free(var_db);
  }
}

lp_variable_db_t* variable_db_new(void) {
  lp_variable_db_t* new = malloc(sizeof(lp_variable_db_t));
  variable_db_construct(new);
  variable_db_attach(new);
  return new;
}

int variable_db_print(const lp_variable_db_t* var_db, FILE* out) {
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

const char* variable_db_get_name(const lp_variable_db_t* var_db, lp_variable_t var) {
  assert(var_db);
  assert(var < var_db->size);
  return var_db->variable_names[var];
}

const lp_variable_db_ops_t lp_variable_db_ops = {
    variable_db_new,
    variable_db_attach,
    variable_db_detach,
    variable_db_new_variable,
    variable_db_add_variable,
    variable_db_print,
    variable_db_get_name
};

#define INITIAL_LIST_SIZE 100
#define INITIAL_MAP_SIZE 100

static
void variable_list_resize(lp_variable_list_t* list, size_t capacity) {
  assert(capacity > list->list_capacity);
  list->list = (lp_variable_t*) realloc(list->list, capacity*sizeof(lp_variable_t));
  list->list_capacity = capacity;
}

static
void variable_map_resize(lp_variable_list_t* list, size_t capacity) {
  assert(capacity > list->var_to_index_map_capacity);
  list->var_to_index_map = (int*) realloc(list->var_to_index_map, capacity*sizeof(int));
  size_t i;
  for (i = list->var_to_index_map_capacity; i < capacity; ++ i) {
    list->var_to_index_map[i] = -1;
  }
  list->var_to_index_map_capacity = capacity;
}

static
void variable_list_construct(lp_variable_list_t* list) {
  // The list
  list->list = 0;
  list->list_size = 0;
  list->list_capacity = 0;
  variable_list_resize(list, INITIAL_LIST_SIZE);
  // Map to indices
  list->var_to_index_map = 0;
  list->var_to_index_map_capacity = 0;
  variable_map_resize(list, INITIAL_MAP_SIZE);
}

static
void variable_list_destruct(lp_variable_list_t* list) {
  free(list->list);
  free(list->var_to_index_map);
}

static
size_t variable_list_size(const lp_variable_list_t* list) {
  return list->list_size;
}

static
int variable_list_index(const lp_variable_list_t* list, lp_variable_t x) {
  if (x >= list->var_to_index_map_capacity) {
    return -1;
  } else {
    return list->var_to_index_map[x];
  }
}

static
void variable_list_copy_into(const lp_variable_list_t* list, lp_variable_t* vars) {
  size_t i;
  for (i = 0; i < list->list_size; ++ i) {
    vars[i] = list->list[i];
  }
}

static
void variable_list_push(lp_variable_list_t* list, lp_variable_t var) {
  if (list->list_size == list->list_capacity) {
    variable_list_resize(list, list->list_capacity*2);
  }
  if (var >= list->var_to_index_map_capacity) {
    variable_map_resize(list, var+1);
  }
  assert(list->var_to_index_map[var] == -1);
  list->var_to_index_map[var] = list->list_size;
  list->list[list->list_size ++] = var;
}

static
void variable_list_pop(lp_variable_list_t* list) {
  assert(list->list_size > 0);
  lp_variable_t var = list->list[-- list->list_size];
  list->var_to_index_map[var] = -1;
}

const lp_variable_list_ops_t lp_variable_list_ops = {
  variable_list_construct,
  variable_list_destruct,
  variable_list_size,
  variable_list_index,
  variable_list_copy_into,
  variable_list_push,
  variable_list_pop
};

