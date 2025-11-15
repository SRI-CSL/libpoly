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
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"
#include "value.h"


/**
 *  Returns the char* representation of the python object pyobj.
 *  If pyobj is NULL, then NULL is returned.
 *  If pyobj is a "string" (PyString, PyBytes, or PyUnicode), then the C string 
 *  associated with the object is returned (using UTF-8).
 *  If pyobj is not a string then we use the str method to obtain
 *  a C string version (again using UTF-8 if necessary).
 *
 */
const char* pythonObject2CharStar(PyObject *pyobj);

/**
 * Check if o is long or int.
 */
int PyLong_or_Int_Check(PyObject* o);

/**
 * Construct the integer from pytong Long or int.
 */
void PyLong_or_Int_to_integer(PyObject* o, const lp_int_ring_t* K, lp_integer_t* x);

/**
 * Get a python integer (or long) from our integer.
 */
PyObject* integer_to_PyInt(const lp_integer_t* x);

/**
 * Get a python long from our integer.
 */
PyObject* integer_to_PyLong(const lp_integer_t* x);

/**
 * Construct the rational from Python float.
 */
void PyFloat_to_dyadic_rational(PyObject* o, lp_dyadic_rational_t* x);

/**
 * Get a python float from our dyadic rational.
 */
PyObject* dyadic_rational_to_PyFloat(const lp_dyadic_rational_t* x);

/**
 * Get a python float from our algebraic number
 */
PyObject* algebraic_number_to_PyFloat(const lp_algebraic_number_t* x);

/**
 * Check if a python number.
 */
int PyLong_or_Int_or_Float_Check(PyObject* o);

/**
 * Construct a value from a python number.
 */
void PyLong_or_Int_or_float_to_value(PyObject* o, lp_value_t* v);
