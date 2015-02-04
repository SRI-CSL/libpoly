/*
 * Coefficient.h
 *
 *  Created on: Nov 8, 2013
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "integer.h"

/** Definition of the python CoefficientRing object */
typedef struct {
  PyObject_HEAD
  lp_int_ring K;
} CoefficientRing;

/** Methods on coefficient rings */
extern PyMethodDef CoefficientRing_methods[];

/** Definition of the CoefficientRing type */
extern PyTypeObject CoefficientRingType;

/** Create a ring object */
PyObject* PyCoefficientRing_create(lp_int_ring K);

/** Check if the object is a polynomial */
#define PyCoefficientRing_CHECK(arg) \
  (Py_TYPE(arg) == &CoefficientRingType)
