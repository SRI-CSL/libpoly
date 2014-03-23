/*
 * utils.h
 *
 *  Created on: Jan 16, 2014
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "integer.h"
#include "rational.h"
#include "dyadic_rational.h"
#include "algebraic_number.h"

/**
 * Check if o is long or int.
 */
int PyLong_or_Int_Check(PyObject* o);

/**
 * Construct the integer from pytong Long or int.
 */
void PyLong_or_Int_to_integer(PyObject* o, int_ring K, integer_t* x);

/**
 * Get a python integer from our integer.
 */
PyObject* integer_to_PyInt(const integer_t* x);

/**
 * Construct the rational from Python float.
 */
void PyFloat_to_dyadic_rational(PyObject* o, dyadic_rational_t* x);

/**
 * Get a python float from our dyadic rational.
 */
PyObject* dyadic_rational_to_PyFloat(const dyadic_rational_t* x);

/**
 * Get a python float from our algebraic number
 */
PyObject* algebraic_number_to_PyFloat(const algebraic_number_t* x);
