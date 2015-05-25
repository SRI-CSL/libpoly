/*
 * FeasibilitySet.h
 *
 *  Created on: May 25, 2015
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "feasibility_set.h"

/** Definition of the python FeasibilitySet object */
typedef struct {
  PyObject_HEAD
  lp_feasibility_set_t* S;
} FeasibilitySet;

/** Methods on intervals */
extern PyMethodDef FeasibilitySet_methods[];

/** Definition of the FeasibilitySet type */
extern PyTypeObject FeasibilitySetType;

/** Create a FeasibilitySet (takes over S) */
PyObject* PyFeasibilitySet_create(lp_feasibility_set_t* S);

/** Check if the object is a polynomial */
#define PyFeasibilitySet_CHECK(arg) \
  (Py_TYPE(arg) == &FeasibilitySetType)
