/*
 * variable.h
 *
 *  Created on: Jan 28, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include <stdio.h>

/** Variables are just ids */
typedef size_t lp_variable_t;

typedef struct {
  /** Reference count */
  size_t ref_count;
  /** Size of the current variable database */
  size_t size;
  /** The capacity of the database */
  size_t capacity;
  /** Names of individual variables */
  char** variable_names;
} lp_variable_db_t;


/** Construct and attach an empty database */
lp_variable_db_t* lp_variable_db_new(void);

/** Attach the variable database */
void lp_variable_db_attach(lp_variable_db_t* var_db);

/** Detach the variable database */
void lp_variable_db_detach(lp_variable_db_t* var_db);

/** Add a new variable to the database and return its id */
lp_variable_t lp_variable_db_new_variable(lp_variable_db_t* var_db, const char* name);

/** Add a new variable to database with a given id*/
void lp_variable_db_add_variable(lp_variable_db_t* var_db, lp_variable_t var, const char* name);

/** Print the variable database */
int lp_variable_db_print(const lp_variable_db_t* var_db, FILE* out);

/** Get the name of the variable */
const char* lp_variable_db_get_name(const lp_variable_db_t* var_db, lp_variable_t var);


/** A list of variable that keeps an index for each variable */
typedef struct {
  /** List of variables in the order */
  lp_variable_t *list;
  /** Size of the list */
  size_t list_size;
  /** Capacity of the list */
  size_t list_capacity;
  /** Map from variables to the index in the list (-1 if not in the order) */
  int* var_to_index_map;
  /** Size of the variable map */
  size_t var_to_index_map_capacity;
} lp_variable_list_t;


/** Construct a new variable order */
void lp_variable_list_construct(lp_variable_list_t* list);

/** Destruct the variable order */
void lp_variable_list_destruct(lp_variable_list_t* list);

/** Get the size of the list */
size_t lp_variable_list_size(const lp_variable_list_t* list);

/** Get the index of the variable in the list, or -1 if not there */
int lp_variable_list_index(const lp_variable_list_t* list, lp_variable_t x);

/** Copy the variables into the given vector */
void lp_variable_list_copy_into(const lp_variable_list_t* list, lp_variable_t* vars);

/** Push a variable to the list */
void lp_variable_list_push(lp_variable_list_t* list, lp_variable_t var);

/** Pop the last variable from the list */
void lp_variable_list_pop(lp_variable_list_t* list);


