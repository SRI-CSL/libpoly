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

#include "integer.h"

/** Definition of the python CoefficientRing object */
typedef struct {
  PyObject_HEAD
  lp_int_ring_t* K;
} CoefficientRing;

/** Methods on coefficient rings */
extern PyMethodDef CoefficientRing_methods[];

/** Definition of the CoefficientRing type */
extern PyTypeObject CoefficientRingType;

/** Create a ring object */
PyObject* PyCoefficientRing_create(lp_int_ring_t* K);

/** Check if the object is a polynomial */
#define PyCoefficientRing_CHECK(arg) \
  (Py_TYPE(arg) == &CoefficientRingType)
