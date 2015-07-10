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

#include "polypyInterval.h"
#include "polypyValue.h"

#include <structmember.h>

static void
Interval_dealloc(Interval* self);

static PyObject*
Interval_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
Interval_init(Interval* self, PyObject* args);

static PyObject*
Interval_str(PyObject* self);

static PyObject*
Interval_pick_value(PyObject* self);

static PyObject*
Interval_contains_value(PyObject* self, PyObject* args);

PyMethodDef Interval_methods[] = {
    {"pick_value", (PyCFunction)Interval_pick_value, METH_NOARGS, "Returns a value from the interval."},
    {"contains", (PyCFunction)Interval_contains_value, METH_VARARGS, "Returns true if the value is in the interval."},
    {NULL}  /* Sentinel */
};

PyTypeObject IntervalType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "polypy.Interval",   /*tp_name*/
    sizeof(Interval), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)Interval_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0                    ,      /*tp_compare*/
    Interval_str,               /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    Interval_str,            /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "Interval objects", /* tp_doc */
    0,                             /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Interval_methods,       /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Interval_init,/* tp_init */
    0,                         /* tp_alloc */
    Interval_new,           /* tp_new */
};

static void
Interval_dealloc(Interval* self)
{
  lp_interval_destruct(&self->I);
  self->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyInterval_create(const lp_interval_t* I) {
  Interval *self;
  self = (Interval*)IntervalType.tp_alloc(&IntervalType, 0);
  if (self != NULL && I) {
    lp_interval_construct_copy(&self->I, I);
  }
  return (PyObject *)self;
}

static PyObject*
Interval_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyInterval_create(0);
}

static int
Interval_init(Interval* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 0) {
    // Defaults to (-inf, +inf)
    lp_interval_construct_full(&self->I);
  } else {
    return -1;
  }
  return 0;
}

static PyObject* Interval_str(PyObject* self) {
  Interval* I = (Interval*) self;
  char* I_str = lp_interval_to_string(&I->I);
  PyObject* str = PyString_FromString(I_str);
  free(I_str);
  return str;
}

static PyObject*
Interval_pick_value(PyObject* self) {
  Interval* I = (Interval*) self;
  lp_value_t v;
  lp_value_construct_none(&v);
  lp_interval_pick_value(&I->I, &v);
  PyObject* result = PyValue_create(&v);
  lp_value_destruct(&v);
  return result;
}

static PyObject*
Interval_contains_value(PyObject* self, PyObject* args) {

  if (!PyTuple_Check(args) || PyTuple_Size(args) != 1) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  PyObject* value_obj = PyTuple_GetItem(args, 0);

  if (!PyValue_CHECK(value_obj)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  lp_interval_t* I = &((Interval*) self)->I;
  lp_value_t* v = &((Value*) value_obj)->v;

  int result = lp_interval_contains(I, v);

  PyObject* result_object = result ? Py_True : Py_False;
  Py_INCREF(result_object);

  return result_object;
}
