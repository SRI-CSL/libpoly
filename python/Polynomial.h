/*
 * Polynomial.h
 *
 *  Created on: Feb 5, 2014
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include "utils.h"
#include "polynomial.h"

/** Definition of the python UPolynomial object */
typedef struct {
  PyObject_HEAD
  lp_polynomial_t* p;
} Polynomial;

/** Methods on univariate polynomials */
extern PyMethodDef Polynomial_methods[];

/** Overloads of arithmetic on polynomial */
extern PyNumberMethods Polynomial_NumberMethods;

/** Definition of the python type */
extern PyTypeObject PolynomialType;

/**
 * Create a new polynomial. It will take over p into the new polynomial hence
 * don't destroy it.
 */
PyObject* Polynomial_create(lp_polynomial_t* p);

/** Get the default context */
const lp_polynomial_context_t* Polynomial_get_default_context(void);

/** Create a polynomial from a variable (new object with refcount 1) */
PyObject* PyPolynomial_FromVariable(PyObject* variable, const lp_polynomial_context_t* ctx);

/** Create a polynomial from int or a long (new object with refcount 1) */
PyObject* PyPolynomial_FromLong_or_Int(PyObject* number, const lp_polynomial_context_t* ctx);

/** Check if the object is a polynomial */
#define PyPolynomial_CHECK(arg) \
  (Py_TYPE(arg) == &PolynomialType)

