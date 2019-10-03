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
