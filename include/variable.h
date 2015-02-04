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

/** Interface to the variable database */
typedef struct {

  /** Construct and attach an empty database */
  lp_variable_db_t* (*new) (void);

  /** Attach the variable database */
  void (*attach) (lp_variable_db_t* var_db);

  /** Detach the variable database */
  void (*detach) (lp_variable_db_t* var_db);

  /** Add a new variable to the database and return its id */
  lp_variable_t (*new_variable) (lp_variable_db_t* var_db, const char* name);

  /** Add a new variable to database with a given id*/
  void (*add_variable) (lp_variable_db_t* var_db, lp_variable_t var, const char* name);

  /** Print the variable database */
  int (*print) (const lp_variable_db_t* var_db, FILE* out);

  /** Get the name of the variable */
  const char* (*get_name) (const lp_variable_db_t* var_db, lp_variable_t var);

} lp_variable_db_ops_t;

/** Implementation of the variable interface */
extern const lp_variable_db_ops_t lp_variable_db_ops;

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

typedef struct {

  /** Construct a new variable order */
  void (*construct) (lp_variable_list_t* list);

  /** Destruct the variable order */
  void (*destruct) (lp_variable_list_t* list);

  /** Get the size of the list */
  size_t (*size) (const lp_variable_list_t* list);

  /** Get the index of the variable in the list, or -1 if not there */
  int (*index) (const lp_variable_list_t* list, lp_variable_t x);

  /** Copy the variables into the given vector */
  void (*copy_into) (const lp_variable_list_t* list, lp_variable_t* vars);

  /** Push a variable to the list */
  void (*push) (lp_variable_list_t* list, lp_variable_t var);

  /** Pop the last variable from the list */
  void (*pop) (lp_variable_list_t* list);

} lp_variable_list_ops_t;

extern const lp_variable_list_ops_t lp_variable_list_ops;
