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

#include "polypyVariableOrder.h"
#include "polypyVariable.h"

#include <structmember.h>
#include <stdlib.h>

/** Default variable database */
static lp_variable_order_t* default_var_order = 0;

lp_variable_order_t* VariableOrder_get_default_order(void) {
  if (!default_var_order) {
    default_var_order = lp_variable_order_new();
  }
  return default_var_order;
}

static PyObject*
VariableOrder_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
VariableOrder_init(VariableOrder* self, PyObject* args);

static void
VariableOrder_dealloc(VariableOrder* self);

static PyObject*
VariableOrder_str(PyObject* self);

static PyObject*
VariableOrder_repr(PyObject* self);

static PyObject*
VariableOrder_push(PyObject* self, PyObject* args);

static PyObject*
VariableOrder_pop(PyObject* self);

static PyObject*
VariableOrder_set(PyObject* self, PyObject* args);

static PyObject*
VariableOrder_cmp(PyObject* self, PyObject* args);

PyMethodDef VariableOrder_methods[] = {
    {"push", (PyCFunction)VariableOrder_push, METH_VARARGS, "Add a variable to the top of the order"},
    {"pop", (PyCFunction)VariableOrder_pop, METH_NOARGS, "Remove the top variable from the order"},
    {"cmp", (PyCFunction)VariableOrder_cmp, METH_VARARGS, "Compare two variables"},
    {"set", (PyCFunction)VariableOrder_set, METH_VARARGS, "Set the order to the given list of variables"},
    {NULL}  /* Sentinel */
};

PyTypeObject VariableOrderType = {
    PyObject_HEAD_INIT(NULL)
    0,                            /*ob_size*/
    "polypy.VariableOrder",       /*tp_name*/
    sizeof(VariableOrder),        /*tp_basicsize*/
    0,                            /*tp_itemsize*/
    (destructor)VariableOrder_dealloc, /*tp_dealloc*/
    0,                            /*tp_print*/
    0,                            /*tp_getattr*/
    0,                            /*tp_setattr*/
    0,                            /*tp_compare*/
    VariableOrder_repr,           /*tp_repr*/
    0,                            /*tp_as_number*/
    0,                            /*tp_as_sequence*/
    0,                            /*tp_as_mapping*/
    0,                            /*tp_hash */
    0,                            /*tp_call*/
    VariableOrder_str,                 /*tp_str*/
    0,                            /*tp_getattro*/
    0,                            /*tp_setattro*/
    0,                            /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "VariableOrder objects",      /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    VariableOrder_methods,             /* tp_methods */
    0,                            /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)VariableOrder_init,      /* tp_init */
    0,                            /* tp_alloc */
    VariableOrder_new,                 /* tp_new */
};

PyObject*
VariableOrder_create(lp_variable_order_t* var_order) {
  VariableOrder *self = (VariableOrder*)VariableOrderType.tp_alloc(&VariableOrderType, 0);
  if (self != NULL) {
    self->var_order = var_order;
  }
  return (PyObject *)self;
}

static PyObject*
VariableOrder_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return VariableOrder_create(0);
}

static PyObject*
VariableOrder_set(PyObject* self, PyObject* args) {

  VariableOrder* order = (VariableOrder*) self;
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* list = PyTuple_GetItem(args, 0);
    if (PyList_Check(list)) {
      int i;
      for (i = 0; i < PyList_Size(list); ++ i) {
        if (!PyVariable_CHECK(PyList_GetItem(list, i))) {
          Py_INCREF(Py_NotImplemented);
          return Py_NotImplemented;
        }
      }
      // Clear the current order
      lp_variable_order_clear(order->var_order);
      // Fill the order
      for (i = 0; i < PyList_Size(list); ++ i) {
        lp_variable_t var = ((Variable*) PyList_GetItem(list, i))->x;
        if (!lp_variable_order_contains(order->var_order, var)) {
          lp_variable_order_push(order->var_order, var);
        }
      }
    } else {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Py_RETURN_NONE;
}

static int
VariableOrder_init(VariableOrder* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* list = PyTuple_GetItem(args, 0);
    if (PyList_Check(list)) {
      int i;
      for (i = 0; i < PyList_Size(list); ++ i) {
        if (!PyVariable_CHECK(PyList_GetItem(list, i))) {
          return -1;
        }
      }
      self->var_order = lp_variable_order_new();
      for (i = 0; i < PyList_Size(list); ++ i) {
        lp_variable_t var = ((Variable*) PyList_GetItem(list, i))->x;
        lp_variable_order_push(self->var_order, var);
      }
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

static void
VariableOrder_dealloc(VariableOrder* self)
{
  if (self->var_order) {
    lp_variable_order_detach(self->var_order);
  }
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
VariableOrder_str(PyObject* self) {
  VariableOrder* var_order = (VariableOrder*) self;
  char* var_order_str = lp_variable_order_to_string(var_order->var_order, Variable_get_default_db());
  PyObject* str = PyString_FromString(var_order_str);
  free(var_order_str);
  return str;
}

static PyObject*
VariableOrder_repr(PyObject* self) {
  VariableOrder* var_order = (VariableOrder*) self;
  char* var_order_str = lp_variable_order_to_string(var_order->var_order, Variable_get_default_db());
  char* var_order_repr = malloc(strlen(var_order_str) + strlen(VariableOrderType.tp_name) + 3);
  sprintf(var_order_repr, "%s(%s)", VariableOrderType.tp_name, var_order_str);
  PyObject* str = PyString_FromString(var_order_repr);
  free(var_order_repr);
  free(var_order_str);
  return str;
}

static PyObject*
VariableOrder_push(PyObject* self, PyObject* args) {
  int error = 0;
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* variable = PyTuple_GetItem(args, 0);
    if (PyVariable_CHECK(variable)) {
      lp_variable_order_t* var_order = ((VariableOrder*) self)->var_order;
      lp_variable_t x = ((Variable*) variable)->x;
      if (!lp_variable_order_contains(var_order, x)) {
        lp_variable_order_push(var_order, x);
      }
    } else {
      error = 1;
    }
  } else {
    error = 1;
  }
  if (error) {
    PyErr_SetString(PyExc_BaseException, "Only variables can be pushed");
  }

  Py_RETURN_NONE;
}


static PyObject*
VariableOrder_cmp(PyObject* self, PyObject* args) {
  int error = 0;
  if (PyTuple_Check(args) && PyTuple_Size(args) == 2) {
    PyObject* x_var = PyTuple_GetItem(args, 0);
    PyObject* y_var = PyTuple_GetItem(args, 1);
    if (PyVariable_CHECK(x_var) && PyVariable_CHECK(y_var)) {
      lp_variable_order_t* var_order = ((VariableOrder*) self)->var_order;
      lp_variable_t x = ((Variable*) x_var)->x;
      lp_variable_t y = ((Variable*) y_var)->x;
      int cmp = lp_variable_order_cmp(var_order, x, y);
      if (cmp < 0) return PyInt_FromLong(-1);
      else if (cmp > 0) return PyInt_FromLong(1);
      else return PyInt_FromLong(0);
    } else {
      error = 1;
    }
  } else {
    error = 1;
  }
  if (error) {
    PyErr_SetString(PyExc_BaseException, "Only variables can be pushed");
  }

  Py_RETURN_NONE;
}


static PyObject*
VariableOrder_pop(PyObject* self) {
  lp_variable_order_t* var_order = ((VariableOrder*) self)->var_order;
  if (lp_variable_order_size(var_order) > 0) {
    lp_variable_order_pop(var_order);
  }
  Py_RETURN_NONE;
}
