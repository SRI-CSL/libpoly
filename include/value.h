/*
 * value.h
 *
 *  Created on: Apr 3, 2014
 *      Author: dejan
 */

#pragma once

#include "poly.h"

#include "integer.h"
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"

#include <stdio.h>

/** Types of values for the assignment */
typedef enum {
  VALUE_NONE,
  VALUE_INTEGER,
  VALUE_DYADIC_RATIONAL,
  VALUE_RATIONAL,
  VALUE_ALGEBRAIC
} value_type_t;

/** A value is a choice of the available types */
typedef union {
  integer_t z;
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
  /** Compare two values. */
  int (*cmp) (const value_t* v1, const value_t* v2);
  /** Void version of the comparison, use with care. */
  int (*cmp_void) (const void* v1, const void* v2);
  /** Print the value */
  int (*print) (const value_t* v, FILE* out);
} value_ops_t;

/** Implementation of the value operations */
extern const value_ops_t value_ops;

