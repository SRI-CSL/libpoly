/*
 * variable_order.h
 *
 *  Created on: Feb 7, 2014
 *      Author: dejan
 */

#pragma once

#include "variable.h"

typedef struct variable_order_struct variable_order_t;
typedef struct variable_order_ops_struct variable_order_ops_t;

typedef struct variable_order_struct {
  variable_order_ops_t* ops;
} variable_order_t;

/**
 * Variable order interface. The order should be total, but can change over
 * time.
 *
 * The interface assumes that the variables are ordered using a score function,
 * that there is an ever-increasing global-timestamp, and that a variable is
 * assigned this time-stamp
 *
 * This can be used to check if a subset of all variables has changed its order.
 * To do so we can keep the maximal time-stamp of all the variables. If the
 * timestamp of any of the variables
 */
typedef struct variable_order_ops_struct {

  /** Construct a new variable order */
  variable_order_t* (*new) (void);

  /** Attach an object to this order (constructor should attach) */
  void (*attach) (variable_order_t* var_order);

  /** Detach an object from this order (frees if refcount = 0) */
  void (*detach) (variable_order_t* var_order);

  /** Compare two variables */
  int (*cmp) (const variable_order_t* var_order, variable_t x, variable_t y);

} variable_order_ops_t;

typedef struct variable_order_simple_struct variable_order_simple_t;
typedef struct variable_order_simple_ops_struct variable_order_simple_ops_t;


/**
 * A simple variable order that orders variable based on a given list, and
 * order the rest of the variables based on their variable id.
 */
typedef struct variable_order_simple_struct {
  /** The operations */
  variable_order_simple_ops_t* ops;
  /** Reference count */
  size_t ref_count;
  /** The actual order */
  variable_list_t list;
} variable_order_simple_t;

typedef struct variable_order_simple_ops_struct {

  variable_order_ops_t variable_order_ops;

  /** Get the size of the order */
  size_t (*size) (const variable_order_simple_t* var_order);

  /** Clear the order */
  void (*clear) (variable_order_simple_t* var_order);

  /** Does the order have an opinion on x */
  int (*contains) (variable_order_simple_t* var_order, variable_t x);

  /** Push a variable to the list */
  void (*push) (variable_order_simple_t* var_order, variable_t var);

  /** Pop the last variable from the list */
  void (*pop) (variable_order_simple_t* var_order);

  /** Print the list of variables */
  int (*print) (const variable_order_simple_t* var_order, const variable_db_t* var_db, FILE* out);

  /** Return a string representation of the order */
  char* (*to_string) (const variable_order_simple_t* var_order, const variable_db_t* var_db);

} variable_order_simple_ops_t;

extern variable_order_simple_ops_t variable_order_simple_ops;
