/*
 * FeasibilitySet.c
 *
 *  Created on: May 25, 2015
 *      Author: dejan
 */

#include "FeasibilitySet.h"
#include "Value.h"

#include <structmember.h>

static void
FeasibilitySet_dealloc(FeasibilitySet* self);

static PyObject*
FeasibilitySet_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
FeasibilitySet_init(FeasibilitySet* self, PyObject* args);

static PyObject*
FeasibilitySet_str(PyObject* self);

static PyObject*
FeasibilitySet_pick_value(PyObject* self);

static PyObject*
FeasibilitySet_intersect(PyObject* self, PyObject* args);

PyMethodDef FeasibilitySet_methods[] = {
    {"pick_value", (PyCFunction)FeasibilitySet_pick_value, METH_NOARGS, "Returns a value from the interval."},
    {"intersect", (PyCFunction)FeasibilitySet_intersect, METH_VARARGS, "Returns the intersection of the interval with another one."},
    {NULL}  /* Sentinel */
};

PyTypeObject FeasibilitySetType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "polypy.FeasibilitySet",   /*tp_name*/
    sizeof(FeasibilitySet), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)FeasibilitySet_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0                    ,      /*tp_compare*/
    FeasibilitySet_str,               /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    FeasibilitySet_str,            /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "FeasibilitySet objects", /* tp_doc */
    0,                             /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    FeasibilitySet_methods,       /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)FeasibilitySet_init,/* tp_init */
    0,                         /* tp_alloc */
    FeasibilitySet_new,           /* tp_new */
};

static void
FeasibilitySet_dealloc(FeasibilitySet* self)
{
  lp_feasibility_set_delete(self->S);
  self->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyFeasibilitySet_create(lp_feasibility_set_t* S) {
  FeasibilitySet *self;
  self = (FeasibilitySet*)FeasibilitySetType.tp_alloc(&FeasibilitySetType, 0);
  if (self != NULL) {
    self->S = S;
  }
  return (PyObject *)self;
}

static PyObject*
FeasibilitySet_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyFeasibilitySet_create(0);
}

static int
FeasibilitySet_init(FeasibilitySet* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 0) {
    // Defaults to (-inf, +inf)
    self->S = lp_feasibility_set_new();
  } else {
    return -1;
  }
  return 0;
}

static PyObject* FeasibilitySet_str(PyObject* self) {
  FeasibilitySet* S = (FeasibilitySet*) self;
  char* I_str = lp_feasibility_set_to_string(S->S);
  PyObject* str = PyString_FromString(I_str);
  free(I_str);
  return str;
}

static PyObject*
FeasibilitySet_pick_value(PyObject* self) {
  FeasibilitySet* S = (FeasibilitySet*) self;
  lp_value_t v;
  lp_value_construct_none(&v);
  lp_feasibility_set_pick_value(S->S, &v);
  PyObject* result = PyValue_create(&v);
  lp_value_destruct(&v);
  return result;
}

static PyObject*
FeasibilitySet_intersect(PyObject* self, PyObject* args) {
  FeasibilitySet* S = (FeasibilitySet*) self;

  if (!PyTuple_Check(args) || PyTuple_Size(args) != 1) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  PyObject* feasibility_set_obj = PyTuple_GetItem(args, 0);
  if (!PyFeasibilitySet_CHECK(feasibility_set_obj)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  // Get the arguments
  lp_feasibility_set_t* S1 = ((FeasibilitySet*) self)->S;
  lp_feasibility_set_t* S2 = ((FeasibilitySet*) feasibility_set_obj)->S;

  // The intersect
  lp_feasibility_set_t* P = lp_feasibility_set_intersect(S1, S2);
  PyObject* P_obj = PyFeasibilitySet_create(P);

  return P_obj;
}
