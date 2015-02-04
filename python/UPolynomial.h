/*
 * UPolynomial.h
 *
 *  Created on: Nov 8, 2013
 *      Author: dejan
 */

#pragma once

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
