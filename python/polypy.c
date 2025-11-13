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

#include "python.h"

#include "polypyInteger.h"
#include "polypyVariable.h"
#include "polypyVariableOrder.h"
#include "polypyUPolynomial.h"
#include "polypyAlgebraicNumber.h"
#include "polypyPolynomial.h"
#include "polypyAssignment.h"
#include "polypyValue.h"
#include "polypyInterval.h"
#include "polypyFeasibilitySet.h"

static PyObject*
Trace_enable(PyObject* self, PyObject* args) {
#ifndef NDEBUG
  const char *tag;
  if (!PyArg_ParseTuple(args, "s", &tag)) {
    return 0;
  }
  lp_trace_enable(tag);
#endif
  Py_RETURN_NONE;
}

static PyObject*
Trace_disable(PyObject* self, PyObject* args) {
#ifndef NDEBUG
  const char *tag;
  if (!PyArg_ParseTuple(args, "s", &tag)) {
    return 0;
  }
  lp_trace_disable(tag);
#endif
  Py_RETURN_NONE;
}

static PyObject*
Stats_print(PyObject* self) {
  lp_stats_print(stdout);
  Py_RETURN_NONE;
}

static PyMethodDef polypy_methods[] = {
    {"trace_enable", (PyCFunction)Trace_enable, METH_VARARGS, "Enables tracing for the given tag"},
    {"trace_disable", (PyCFunction)Trace_disable, METH_VARARGS, "Disables tracing for the given tag"},
    {"stats_print", (PyCFunction)Stats_print, METH_NOARGS, "Prints the statistics"},
    {NULL}  /* Sentinel */
};


static struct PyModuleDef polypymodule = {
  PyModuleDef_HEAD_INIT,
  "polypy",
  NULL,
  0, // sizeof polypy struct
  polypy_methods,
  0,
  0,
  0
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
PyInit_polypy(void)
{
  PyObject* m;

  if (PyType_Ready(&CoefficientRingType) < 0)
    return NULL;
  if (PyType_Ready(&VariableType) < 0)
    return NULL;
  if (PyType_Ready(&VariableOrderType) < 0)
    return NULL;
  if (PyType_Ready(&AssignmentType) < 0)
    return NULL;
  if (PyType_Ready(&PolynomialType) < 0)
    return NULL;
  if (PyType_Ready(&UPolynomialType) < 0)
    return NULL;
  if (PyType_Ready(&AlgebraicNumberType) < 0)
    return NULL;
  if (PyType_Ready(&ValueType) < 0)
    return NULL;
  if (PyType_Ready(&IntervalType) < 0)
    return NULL;
  if (PyType_Ready(&FeasibilitySetType) < 0)
    return NULL;

  m = PyModule_Create(&polypymodule);
  
  // Initialize the library
  lp_set_output_language(LP_OUTPUT_PYTHON);

  Py_INCREF(&CoefficientRingType);
  PyModule_AddObject(m, "CoefficientRing", (PyObject*)&CoefficientRingType);

  PyObject* PyZ = PyCoefficientRing_create(lp_Z);
  Py_INCREF(PyZ);
  PyModule_AddObject(m, "Z", PyZ);

  Py_INCREF(&VariableType);
  PyModule_AddObject(m, "Variable", (PyObject*)&VariableType);

  Py_INCREF(&VariableOrderType);
  PyModule_AddObject(m, "VariableOrder", (PyObject*)&VariableOrderType);

  PyObject* variable_order = VariableOrder_create(VariableOrder_get_default_order());
  Py_INCREF(variable_order);
  PyModule_AddObject(m, "variable_order", variable_order);

  Py_INCREF(&PolynomialType);
  PyModule_AddObject(m, "Polynomial", (PyObject*)&PolynomialType);

  // Sign conditions
  PyObject* Py_SGN_LT_0 = PyLong_FromLong(LP_SGN_LT_0);
  PyObject* Py_SGN_LE_0 = PyLong_FromLong(LP_SGN_LE_0);
  PyObject* Py_SGN_EQ_0 = PyLong_FromLong(LP_SGN_EQ_0);
  PyObject* Py_SGN_NE_0 = PyLong_FromLong(LP_SGN_NE_0);
  PyObject* Py_SGN_GT_0 = PyLong_FromLong(LP_SGN_GT_0);
  PyObject* Py_SGN_GE_0 = PyLong_FromLong(LP_SGN_GE_0);
  PyModule_AddObject(m, "SGN_LT_0", Py_SGN_LT_0);
  PyModule_AddObject(m, "SGN_LE_0", Py_SGN_LE_0);
  PyModule_AddObject(m, "SGN_EQ_0", Py_SGN_EQ_0);
  PyModule_AddObject(m, "SGN_NE_0", Py_SGN_NE_0);
  PyModule_AddObject(m, "SGN_GT_0", Py_SGN_GT_0);
  PyModule_AddObject(m, "SGN_GE_0", Py_SGN_GE_0);

  Py_INCREF(&AssignmentType);
  PyModule_AddObject(m, "Assignment", (PyObject*)&AssignmentType);

  Py_INCREF(&UPolynomialType);
  PyModule_AddObject(m, "UPolynomial", (PyObject*)&UPolynomialType);

  int x_coeff[2] = { 0, 1 };
  lp_upolynomial_t* x_poly = lp_upolynomial_construct_from_int(lp_Z, 1, x_coeff);
  PyObject* x = PyUPolynomial_create(x_poly);
  Py_INCREF(x);
  PyModule_AddObject(m, "x", x);

  Py_INCREF(&AlgebraicNumberType);
  PyModule_AddObject(m, "AlgebraicNumber", (PyObject*)&AlgebraicNumberType);

  Py_INCREF(&ValueType);
  PyModule_AddObject(m, "Value", (PyObject*)&ValueType);

  // Sign conditions
  lp_value_t value_inf_pos, value_inf_neg;
  lp_value_construct(&value_inf_pos, LP_VALUE_PLUS_INFINITY, NULL);
  lp_value_construct(&value_inf_neg, LP_VALUE_MINUS_INFINITY, NULL);
  PyObject* inf_pos = PyValue_create(&value_inf_pos);
  PyObject* inf_neg = PyValue_create(&value_inf_neg);
  PyModule_AddObject(m, "INFINITY_POS", inf_pos);
  PyModule_AddObject(m, "INFINITY_NEG", inf_neg);
  Py_INCREF(inf_pos);
  Py_INCREF(inf_neg);
  lp_value_destruct(&value_inf_pos);
  lp_value_destruct(&value_inf_neg);

  Py_INCREF(&IntervalType);
  PyModule_AddObject(m, "Interval", (PyObject*)&IntervalType);

  Py_INCREF(&FeasibilitySetType);
  PyModule_AddObject(m, "FeasibilitySet", (PyObject*)&FeasibilitySetType);

  return m;
}
