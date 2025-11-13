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

#include "polypyFeasibilitySet.h"
#include "polypyValue.h"

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

static PyObject*
FeasibilitySet_contains_value(PyObject* self, PyObject* args);

static PyObject*
FeasibilitySet_contains_int(PyObject* self);

PyMethodDef FeasibilitySet_methods[] = {
    {"pick_value", (PyCFunction)FeasibilitySet_pick_value, METH_NOARGS, "Returns a value from the interval."},
    {"intersect", (PyCFunction)FeasibilitySet_intersect, METH_VARARGS, "Returns the intersection of the set with another one."},
    {"contains", (PyCFunction)FeasibilitySet_contains_value, METH_VARARGS, "Returns true if the value is in the feasibility set."},
    {"contains_int", (PyCFunction)FeasibilitySet_contains_int, METH_NOARGS, "Returns true if the set contains an integer value."},
    {NULL}  /* Sentinel */
};

PyTypeObject FeasibilitySetType = {
    {PyObject_HEAD_INIT(NULL)},         // PyObject_VAR_HEAD
    "polypy.FeasibilitySet",            // const char *tp_name;
    sizeof(FeasibilitySet),             // Py_ssize_t tp_basicsize;
    0,                                  // Py_ssize_t tp_itemsize;
    (destructor)FeasibilitySet_dealloc, // destructor tp_dealloc;
    0,                                  // printfunc tp_print;
    NULL,                               // getattrfunc tp_getattr;
    NULL,                               // setattrfunc tp_setattr;
    NULL,                               // PyAsyncMethods *tp_as_async;
    FeasibilitySet_str,                 // reprfunc tp_repr;
    NULL,                               // PyNumberMethods *tp_as_number;
    NULL,                               // PySequenceMethods *tp_as_sequence;
    NULL,                               // PyMappingMethods *tp_as_mapping;
    NULL,                               // hashfunc tp_hash;
    NULL,                               // ternaryfunc tp_call;
    FeasibilitySet_str,                 // reprfunc tp_str;
    NULL,                               // getattrofunc tp_getattro;
    NULL,                               // setattrofunc tp_setattro;
    NULL,                               // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,                 // unsigned long tp_flags;
    "FeasibilitySet objects",           // const char *tp_doc;
    NULL,                               // traverseproc tp_traverse;
    NULL,                               // inquiry tp_clear;
    NULL,                               // richcmpfunc tp_richcompare;
    0,                                  // Py_ssize_t tp_weaklistoffset;
    NULL,                               // getiterfunc tp_iter;
    NULL,                               // iternextfunc tp_iternext;
    FeasibilitySet_methods,             // struct PyMethodDef *tp_methods;
    NULL,                               // struct PyMemberDef *tp_members;
    NULL,                               // struct PyGetSetDef *tp_getset;
    NULL,                               // struct _typeobject *tp_base;
    NULL,                               // PyObject *tp_dict;
    NULL,                               // descrgetfunc tp_descr_get;
    NULL,                               // descrsetfunc tp_descr_set;
    0,                                  // Py_ssize_t tp_dictoffset;
    (initproc)FeasibilitySet_init,      // initproc tp_init;
    NULL,                               // allocfunc tp_alloc;
    FeasibilitySet_new,                 // newfunc tp_new;
    NULL,                               // freefunc tp_free;
    NULL,                               // inquiry tp_is_gc;
    NULL,                               // PyObject *tp_bases;
    NULL,                               // PyObject *tp_mro;
    NULL,                               // PyObject *tp_cache;
    NULL,                               // PyObject *tp_subclasses;
    NULL,                               // PyObject *tp_weaklist;
    NULL,                               // destructor tp_del;
    0,                                  // unsigned int tp_version_tag;
    NULL,                               // destructor tp_finalize;
};

static void
FeasibilitySet_dealloc(FeasibilitySet* self) {
  lp_feasibility_set_delete(self->S);
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
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
  return PyFeasibilitySet_create(NULL);
}

static int
FeasibilitySet_init(FeasibilitySet* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 0) {
    // Defaults to (-inf, +inf)
    self->S = lp_feasibility_set_new_full();
  } else {
    return -1;
  }
  return 0;
}

static PyObject* FeasibilitySet_str(PyObject* self) {
  FeasibilitySet* S = (FeasibilitySet*) self;
  char* I_str = lp_feasibility_set_to_string(S->S);
  PyObject* str = PyUnicode_FromString(I_str);
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

  // The intersection
  lp_feasibility_set_t* P = lp_feasibility_set_intersect(S1, S2);
  PyObject* P_obj = PyFeasibilitySet_create(P);

  return P_obj;
}

static PyObject*
FeasibilitySet_contains_value(PyObject* self, PyObject* args) {

  if (!PyTuple_Check(args) || PyTuple_Size(args) != 1) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  PyObject* value_obj = PyTuple_GetItem(args, 0);

  if (!PyValue_CHECK(value_obj)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  lp_feasibility_set_t* S = ((FeasibilitySet*) self)->S;
  lp_value_t* v = &((Value*) value_obj)->v;

  int result = lp_feasibility_set_contains(S, v);

  PyObject* result_object = result ? Py_True : Py_False;
  Py_INCREF(result_object);
  return result_object;
}

static PyObject*
FeasibilitySet_contains_int(PyObject* self) {
  lp_feasibility_set_t* S = ((FeasibilitySet*) self)->S;
  PyObject* result_object = lp_feasibility_set_contains_int(S) ? Py_True : Py_False;
  Py_INCREF(result_object);
  return result_object;
}
