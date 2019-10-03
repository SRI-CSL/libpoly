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

#define PY_SSIZE_T_CLEAN
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
