#include "Integer.h"
#include "Variable.h"
#include "VariableOrder.h"
#include "UPolynomial.h"
#include "AlgebraicNumber.h"
#include "Polynomial.h"
#include "Assignment.h"
#include "Value.h"
#include "Interval.h"

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

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initpolypy(void)
{
  PyObject* m;

  if (PyType_Ready(&CoefficientRingType) < 0)
    return;
  if (PyType_Ready(&VariableType) < 0)
    return;
  if (PyType_Ready(&VariableOrderType) < 0)
    return;
  if (PyType_Ready(&AssignmentType) < 0)
    return;
  if (PyType_Ready(&PolynomialType) < 0)
    return;
  if (PyType_Ready(&UPolynomialType) < 0)
    return;
  if (PyType_Ready(&AlgebraicNumberType) < 0)
    return;
  if (PyType_Ready(&ValueType) < 0)
    return;
  if (PyType_Ready(&IntervalType) < 0)
    return;

  m = Py_InitModule3("polypy", polypy_methods, "PolyPy Libarary.");

  // Initialize the library
  lp_init();
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
  PyObject* Py_SGN_LT_0 = PyInt_FromLong(LP_SGN_LT_0);
  PyObject* Py_SGN_LE_0 = PyInt_FromLong(LP_SGN_LE_0);
  PyObject* Py_SGN_EQ_0 = PyInt_FromLong(LP_SGN_EQ_0);
  PyObject* Py_SGN_NE_0 = PyInt_FromLong(LP_SGN_NE_0);
  PyObject* Py_SGN_GT_0 = PyInt_FromLong(LP_SGN_GT_0);
  PyObject* Py_SGN_GE_0 = PyInt_FromLong(LP_SGN_GE_0);
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
  PyModule_AddObject(m, "Value", (PyObject*)&AlgebraicNumberType);

  Py_INCREF(&IntervalType);
  PyModule_AddObject(m, "Interval", (PyObject*)&AlgebraicNumberType);
}

