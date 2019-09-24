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

#include "polypyInteger.h"

#include <structmember.h>

static void
CoefficientRing_dealloc(CoefficientRing* self);

static PyObject*
CoefficientRing_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
CoefficientRing_init(CoefficientRing* self, PyObject* args);

static int
CoefficientRing_cmp(PyObject* self, PyObject* args);

static PyObject*
CoefficientRing_modulus(PyObject* self);

static PyObject*
CoefficientRing_str(PyObject* self);

PyMethodDef CoefficientRing_methods[] = {
    {"modulus", (PyCFunction)CoefficientRing_modulus, METH_NOARGS, "Returns the degree of the polynomial"},
    {NULL}  /* Sentinel */
};

PyTypeObject CoefficientRingType = {
    PyObject_HEAD_INIT(NULL)
    "polypy.CoefficientRing", //const char *tp_name; /* For printing, in format "<module>.<name>" */
    sizeof(CoefficientRing), //Py_ssize_t tp_basicsize;
    0, //Py_ssize_t tp_itemsize; /* For allocation */
    (destructor)CoefficientRing_dealloc, //destructor tp_dealloc;
    0, // printfunc tp_print;
    0, // getattrfunc tp_getattr;
    0, // setattrfunc tp_setattr;
    0, //CoefficientRing_cmp, //PyAsyncMethods *tp_as_async; /* formerly known as tp_compare (Python 2) or tp_reserved (Python 3) */
    0, // reprfunc tp_repr;
    0, // PyNumberMethods *tp_as_number;
    0, // PySequenceMethods *tp_as_sequence;
    0, // PyMappingMethods *tp_as_mapping;
    0, // hashfunc tp_hash;
    0, // ternaryfunc tp_call;
    CoefficientRing_str, //reprfunc tp_str;
    0, // getattrofunc tp_getattro;
    0, // setattrofunc tp_setattro;
    0, // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT, //unsigned long tp_flags;
    "Coefficient ring objects", //const char *tp_doc; /* Documentation string */
    0, // traverseproc tp_traverse;
    0, // inquiry tp_clear;
    0, // richcmpfunc tp_richcompare;
    0, // Py_ssize_t tp_weaklistoffset;
    0, // getiterfunc tp_iter;
    0, // iternextfunc tp_iternext;
    CoefficientRing_methods, //struct PyMethodDef *tp_methods;
    0, //struct PyMemberDef *tp_members;
    0, //struct PyGetSetDef *tp_getset;
    0, //struct _typeobject *tp_base;
    0, //PyObject *tp_dict;
    0, //descrgetfunc tp_descr_get;
    0, //descrsetfunc tp_descr_set;
    0, //Py_ssize_t tp_dictoffset;
    (initproc)CoefficientRing_init, //initproc tp_init;
    0, //allocfunc tp_alloc;
    CoefficientRing_new, //newfunc tp_new;
    0, //freefunc tp_free; /* Low-level free-memory routine */
    0, //inquiry tp_is_gc; /* For PyObject_IS_GC */
    0, //PyObject *tp_bases;
    0, //PyObject *tp_mro; /* method resolution order */
    0, //PyObject *tp_cache;
    0, //PyObject *tp_subclasses;
    0, //PyObject *tp_weaklist;
    0, //destructor tp_del;
    0, //unsigned int tp_version_tag;
    0, //destructor tp_finalize;
};

static void
CoefficientRing_dealloc(CoefficientRing* self)
{
  if (self->K) lp_int_ring_detach(self->K);
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyCoefficientRing_create(lp_int_ring_t* K) {
  CoefficientRing *self;
  self = (CoefficientRing*)CoefficientRingType.tp_alloc(&CoefficientRingType, 0);
  if (self != NULL) {
    self->K = K;
  }
  return (PyObject *)self;
}

static PyObject*
CoefficientRing_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyCoefficientRing_create(0);
}


static int
CoefficientRing_init(CoefficientRing* self, PyObject* args)
{
    if (PyTuple_Check(args)) {
      if (PyTuple_Size(args) == 0) {
        // Defaults to Z
        self->K = lp_Z;
      } else if (PyTuple_Size(args) == 1) {
        // Get the list of coefficients
        PyObject* modulus = PyTuple_GetItem(args, 0);
        if (!PyLong_Check(modulus)) {
          if (!PyLong_Check(modulus)) {
            return -1;
          } else {
            long M_int = PyLong_AsLong(modulus);
            lp_integer_t M;
            lp_integer_construct_from_int(lp_Z, &M, M_int);
            int is_prime = lp_integer_is_prime(&M);
            self->K = lp_int_ring_create(&M, is_prime);
            lp_integer_destruct(&M);
          }
        } else {
          int overflow = 0;
          long M_int = PyLong_AsLongAndOverflow(modulus, &overflow);
          if (overflow) {
            PyObject* M_str = PyObject_Str(modulus);
            char* M_cstr = PyBytes_AS_STRING(M_str); //IAM: this is a mess; will fix when I find a decent guide
            lp_integer_t M;
            lp_integer_construct_from_string(lp_Z, &M, M_cstr, 10);
            int is_prime = lp_integer_is_prime(&M);
            self->K = lp_int_ring_create(&M, is_prime);
            Py_DECREF(M_str);
          } else if (M_int > 0) {
            lp_integer_t M;
            lp_integer_construct_from_int(lp_Z, &M, M_int);
            int is_prime = lp_integer_is_prime(&M);
            self->K = lp_int_ring_create(&M, is_prime);
          } else {
            return -1;
          }
        }
      } else {
        return -1;
      }
    } else {
      return -1;
    }

    return 0;
}

static int
CoefficientRing_cmp(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyCoefficientRing_CHECK(self) || !PyCoefficientRing_CHECK(other)) {
    // should return -1 and set an exception condition when an error occurred
    return -1;
  }
  // Get arguments
  CoefficientRing* K1 = (CoefficientRing*) self;
  CoefficientRing* K2 = (CoefficientRing*) other;
  // Are they equal
  if (K1->K == K2->K) {
    return 0;
  }
  // Is one of them Z
  if (K1->K == lp_Z) {
    return 1;
  }
  if (K2->K == lp_Z) {
    return -1;
  }
  // Compare
  return lp_integer_cmp(lp_Z, &K1->K->M, &K2->K->M);
}

static PyObject*
CoefficientRing_modulus(PyObject* self) {
  CoefficientRing* K = (CoefficientRing*) self;
  if (K && K->K) {
    char* K_str = lp_integer_to_string(&K->K->M);
    char* p = 0;
    PyObject* M = PyLong_FromString(K_str, &p, 10);
    free(K_str);
    return M;
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject* CoefficientRing_str(PyObject* self) {
  CoefficientRing* K = (CoefficientRing*) self;
  if (K) {
    if (K->K) {
      char* K_str = lp_int_ring_to_string(K->K);
      PyObject* str = PyUnicode_FromString(K_str);
      free(K_str);
      return str;
    } else {
      return PyUnicode_FromString("Z");
    }
  } else {
    Py_RETURN_NONE;
  }
}
