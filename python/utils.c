/*
 * utils.c
 *
 *  Created on: Jan 16, 2014
 *      Author: dejan
 */

#include "utils.h"

int PyLong_or_Int_Check(PyObject* o) {
  if (PyInt_Check(o)) {
    return 1;
  }
  if (PyLong_Check(o)) {
    return 1;
  }
  return 0;
}

void PyLong_or_Int_to_integer(PyObject* o, lp_int_ring_t* K, lp_integer_t* c) {
  if (PyInt_Check(o)) {
    long c_long  = PyInt_AsLong(o);
    lp_integer_construct_from_int(K, c, c_long);
  }
  if (PyLong_Check(o)) {
    PyObject* o_str = PyObject_Str(o);
    char* o_cstr = PyString_AsString(o_str);
    lp_integer_construct_from_string(K, c, o_cstr, 10);
    Py_DECREF(o_str);
  }
}

PyObject* integer_to_PyInt(const lp_integer_t* x) {
  char* str = lp_integer_to_string(x);
  char* str_p = 0;
  PyObject* result = PyInt_FromString(str, &str_p, 10);
  free(str);
  return result;
}

void PyFloat_to_dyadic_rational(PyObject* o, lp_dyadic_rational_t* x) {
  double o_double = PyFloat_AsDouble(o);
  lp_dyadic_rational_construct_from_double(x, o_double);
}

PyObject* dyadic_rational_to_PyFloat(const lp_dyadic_rational_t* x) {
  double x_double = lp_dyadic_rational_to_double(x);
  return PyFloat_FromDouble(x_double);
}


PyObject* algebraic_number_to_PyFloat(const lp_algebraic_number_t* x) {
  double x_double = lp_algebraic_number_to_double(x);
  return PyFloat_FromDouble(x_double);
}
