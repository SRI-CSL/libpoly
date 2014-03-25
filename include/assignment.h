/*
 * assignment.h
 *
 *  Created on: Jan 28, 2014
 *      Author: dejan
 */

#pragma once

#include "variable.h"
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"

#include "polynomial.h"

#include <stdio.h>

/** Types of values for the assignment */
typedef enum {
  VALUE_NONE,
  VALUE_RATIONAL,
  VALUE_DYADIC_RATIONAL,
  VALUE_ALGEBRAIC
} value_type_t;

/** A value is a choice of the available types */
typedef union {
  rational_t q;
  dyadic_rational_t dy_q;
  algebraic_number_t a;
} value_union_t;

/** A value is a tagged union of available type */
typedef struct {
  value_type_t type;
  value_union_t value;
} value_t;

typedef struct {
  /** Construct a value */
  void (*construct) (value_t* v, value_type_t type, const void* data);
  /** Construct a copy of the given value */
  void (*construct_copy) (value_t* v, const value_t* from);
  /** Destruct the value */
  void (*destruct) (value_t* v);
  /** Get the approximate value */
  void (*approximate) (const value_t* v, interval_t* approx);
  /** Print the value */
  int (*print) (const value_t* v, FILE* out);
} value_ops_t;

/** Implementation of the value operations */
extern const value_ops_t value_ops;

typedef struct {
  /** Size of the map */
  size_t size;
  /** The values */
  value_t* values;
  /** The variable database */
  const variable_db_t* var_db;
} assignment_t;

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
