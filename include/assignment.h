/*
 * assignment.h
 *
 *  Created on: Jan 28, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "variable.h"
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"

#include "polynomial.h"
#include "value.h"

#include <stdio.h>

struct assignment_struct {
  /** Size of the map */
  size_t size;
  /** The values */
  value_t* values;
  /** The variable database */
  const variable_db_t* var_db;
};

typedef struct {

  /** Construct an empty assignment */
  void (*construct) (assignment_t* m, const variable_db_t* var_db);

  /** Construct an empty assignment */
  assignment_t* (*new) (const variable_db_t* var_db);

  /** Destruct the assignment */
  void (*destruct) (assignment_t* m);

  /** Print the model */
  int (*print) (const assignment_t* m, FILE* out);

  /** Get the string representation of the model */
  char* (*to_string) (const assignment_t* m);

  /**
   * Set the value of a variable (value is copied over). If value is 0 (pointer)
   * the value is unset.
   */
  void (*set_value) (assignment_t* m, variable_t x, const value_t* value);

  /** Get the value of a variable */
  const value_t* (*get_value) (const assignment_t* m, variable_t x);

  /** Get an approximate value of the variable */
  void (*get_value_approx) (const assignment_t* m, variable_t x, interval_t* approx);

  /** Get the sign of the polynomial in the model */
  int (*sgn) (const assignment_t* m, const polynomial_t* A);

} assignment_ops_t;

/** Implementation of the assignment interface */
extern const assignment_ops_t assignment_ops;
