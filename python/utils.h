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
#include "value.h"

/**
 * Check if o is long or int.
 */
int PyLong_or_Int_Check(PyObject* o);

/**
 * Construct the integer from pytong Long or int.
 */
void PyLong_or_Int_to_integer(PyObject* o, lp_int_ring_t* K, lp_integer_t* x);

/**
 * Get a python integer from our integer.
 */
PyObject* integer_to_PyInt(const lp_integer_t* x);

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
