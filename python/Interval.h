/*
 * Interval.h
 *
 *  Created on: May 19, 2015
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "interval.h"

/** Definition of the python Interval object */
typedef struct {
  PyObject_HEAD
  lp_interval_t I;
} Interval;

/** Methods on intervals */
extern PyMethodDef Interval_methods[];

/** Definition of the Interval type */
extern PyTypeObject IntervalType;

/** Create an Interval, (makes a copy of I) */
PyObject* PyInterval_create(const lp_interval_t* I);

/** Check if the object is a polynomial */
#define PyInterval_CHECK(arg) \
  (Py_TYPE(arg) == &IntervalType)
