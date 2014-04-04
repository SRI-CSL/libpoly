#include "Integer.h"
#include "Variable.h"
#include "VariableOrder.h"
#include "UPolynomial.h"
#include "AlgebraicNumber.h"
#include "Polynomial.h"
#include "Assignment.h"

static PyObject*
Trace_enable(PyObject* self, PyObject* args) {
#ifndef NDEBUG
  const char *tag;
  if (!PyArg_ParseTuple(args, "s", &tag)) {
    return 0;
  }
  poly_ops.trace_enable(tag);
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
  poly_ops.trace_disable(tag);
#endif
  Py_RETURN_NONE;
}

static PyObject*
Stats_print(PyObject* self) {
  poly_ops.stats_print(stdout);
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

  m = Py_InitModule3("polypy", polypy_methods, "PolyPy Libarary.");

  // Set ** as the power symbol
  polynomial_ops.set_power_symbol("**");

  Py_INCREF(&CoefficientRingType);
  PyModule_AddObject(m, "CoefficientRing", (PyObject*)&CoefficientRingType);

  PyObject* PyZ = PyCoefficientRing_create(Z);
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

  Py_INCREF(&AssignmentType);
  PyModule_AddObject(m, "Assignment", (PyObject*)&AssignmentType);

  Py_INCREF(&UPolynomialType);
  PyModule_AddObject(m, "UPolynomial", (PyObject*)&UPolynomialType);

  int x_coeff[2] = { 0, 1 };
  upolynomial_t* x_poly = upolynomial_ops.construct_from_int(Z, 1, x_coeff);
  PyObject* x = PyUPolynomial_create(x_poly);
  Py_INCREF(x);
  PyModule_AddObject(m, "x", x);

  Py_INCREF(&AlgebraicNumberType);
  PyModule_AddObject(m, "AlgebraicNumber", (PyObject*)&AlgebraicNumberType);

}

