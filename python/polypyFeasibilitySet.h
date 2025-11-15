/**
 * Copyright 2015, SRI International.
 *
 * This file is part of LibPoly.
 *
 * LibPoly is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LibPoly is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibPoly.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "python.h"

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
