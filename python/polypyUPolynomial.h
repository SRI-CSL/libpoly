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

#include "utils.h"
#include "upolynomial.h"
#include "upolynomial_factors.h"

/** Definition of the python UPolynomial object */
typedef struct {
  PyObject_HEAD
  lp_upolynomial_t* p;
} UPolynomialObject;

/** Methods on univariate polynomials */
extern PyMethodDef UPolynomial_methods[];

/** Overloads of arithmetic on polynomial */
extern PyNumberMethods UPolynomial_NumberMethods;

/** Definition of the python type */
extern PyTypeObject UPolynomialType;

/** Create a new polynomial */
PyObject* PyUPolynomial_create(lp_upolynomial_t* p);

/** Check if the object is a polynomial */
#define PyUPolynomial_CHECK(arg) \
  (Py_TYPE(arg) == &UPolynomialType && ((UPolynomialObject*)arg)->p)
