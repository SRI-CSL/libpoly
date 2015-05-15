/*
 * variable_list.h
 *
 *  Created on: Feb 5, 2015
 *      Author: dejan
 */

#pragma once

#include "poly.h"

/**
 * A list of variable that keeps an index for each variable. We don't expect
 * too many variables in the database, so we keep a map as an array.
 */
struct lp_variable_list_struct {
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
};

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

/** Get the last variable from the list */
lp_variable_t lp_variable_list_top(const lp_variable_list_t* list);
