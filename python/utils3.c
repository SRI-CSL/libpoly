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

#include "utils.h"

int PyLong_or_Int_Check(PyObject* o) {
  if (PyLong_Check(o)) {
    return 1;
  }
  return 0;
}

void PyLong_or_Int_to_integer(PyObject* o, lp_int_ring_t* K, lp_integer_t* c) {
  if (PyLong_Check(o)) {
    long c_long  = PyLong_AsLong(o);
    lp_integer_construct_from_int(K, c, c_long);
  }
  if (PyLong_Check(o)) {
    PyObject* o_str = PyObject_Str(o);
    char* o_cstr = PyBytes_AS_STRING(o_str);
    lp_integer_construct_from_string(K, c, o_cstr, 10);
    Py_DECREF(o_str);
  }
}

//IAM: unused in the python3 world
PyObject* integer_to_PyInt(const lp_integer_t* x) {
  char* str = lp_integer_to_string(x);
  char* str_p = 0;
  PyObject* result = PyLong_FromString(str, &str_p, 10);
  free(str);
  return result;
}

PyObject* integer_to_PyLong(const lp_integer_t* x) {
  char* str = lp_integer_to_string(x);
  char* str_p = 0;
  PyObject* result = PyLong_FromString(str, &str_p, 10);
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

int PyLong_or_Int_or_Float_Check(PyObject* o) {
  if (PyLong_Check(o)) {
    return 1;
  }
  if (PyLong_Check(o)) {
    return 1;
  }
  if (PyFloat_Check(o)) {
    return 1;
  }
  return 0;
}

void PyLong_or_Int_or_float_to_value(PyObject* o, lp_value_t* v) {
  if (PyLong_or_Int_Check(o)) {
    lp_integer_t v_int;
    PyLong_or_Int_to_integer(o, lp_Z, &v_int);
    lp_value_construct(v, LP_VALUE_INTEGER, &v_int);
    lp_integer_destruct(&v_int);
  } else if (PyFloat_Check(o)) {
    lp_dyadic_rational_t v_dy_q;
    PyFloat_to_dyadic_rational(o, &v_dy_q);
    lp_value_construct(v, LP_VALUE_INTEGER, &v_dy_q);
    lp_dyadic_rational_destruct(&v_dy_q);
  } else {
    lp_value_construct(v, LP_VALUE_NONE, 0);
  }
}
