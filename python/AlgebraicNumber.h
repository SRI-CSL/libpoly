/*
 * AlgebraicNumber.h
 *
 *  Created on: Jan 26, 2014
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "algebraic_number.h"

/** Definition of the python AlgebraicNumber object */
typedef struct {
  PyObject_HEAD
  algebraic_number_t a;
} AlgebraicNumber;

/** Methods on coefficient rings */
extern PyMethodDef AlgebraicNumber_methods[];

/** Definition of the CoefficientRing type */
extern PyTypeObject AlgebraicNumberType;

/** Create an algebraic object (makes a copy of a) */
PyObject* PyAlgebraicNumber_create(const algebraic_number_t* a);

/** Check if the object is a polynomial */
#define PyAlgebraicNumber_CHECK(arg) \
  (Py_TYPE(arg) == &AlgebraicNumberType)
