/*
 * Value.h
 *
 *  Created on: Jan 26, 2014
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "value.h"

/** Definition of the python Value object */
typedef struct {
  PyObject_HEAD
  lp_value_t v;
} Value;

/** Methods on coefficient rings */
extern PyMethodDef Value_methods[];

/** Definition of the CoefficientRing type */
extern PyTypeObject ValueType;

/** Create an algebraic object (makes a copy of v) */
PyObject* PyValue_create(const lp_value_t* v);

/** Check if the object is a polynomial */
#define PyValue_CHECK(arg) \
  (Py_TYPE(arg) == &ValueType)
