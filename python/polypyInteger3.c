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
#include "utils.h"

#include <structmember.h>

static void
CoefficientRing_dealloc(CoefficientRing* self);

static PyObject*
CoefficientRing_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
CoefficientRing_init(CoefficientRing* self, PyObject* args);

static PyObject*
CoefficientRing_modulus(PyObject* self);

static PyObject*
CoefficientRing_str(PyObject* self);

static PyObject *
CoefficientRing_richcompare(PyObject *self, PyObject *other, int op);

PyMethodDef CoefficientRing_methods[] = {
    {"modulus", (PyCFunction)CoefficientRing_modulus, METH_NOARGS, "Returns the degree of the polynomial"},
    {NULL}  /* Sentinel */
};

PyTypeObject CoefficientRingType = {
    {PyObject_HEAD_INIT(NULL)},          // PyObject_VAR_HEAD
    "polypy.CoefficientRing",            // const char *tp_name;
    sizeof(CoefficientRing),             // Py_ssize_t tp_basicsize;
    0,                                   // Py_ssize_t tp_itemsize;
    (destructor)CoefficientRing_dealloc, // destructor tp_dealloc;
    0,                                   // printfunc tp_print;
    0,                                   // getattrfunc tp_getattr;
    0,                                   // setattrfunc tp_setattr;
    0,                                   // PyAsyncMethods *tp_as_async;
    0,                                   // reprfunc tp_repr;
    0,                                   // PyNumberMethods *tp_as_number;
    0,                                   // PySequenceMethods *tp_as_sequence;
    0,                                   // PyMappingMethods *tp_as_mapping;
    0,                                   // hashfunc tp_hash;
    0,                                   // ternaryfunc tp_call;
    CoefficientRing_str,                 // reprfunc tp_str;
    0,                                   // getattrofunc tp_getattro;
    0,                                   // setattrofunc tp_setattro;
    0,                                   // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,                  // unsigned long tp_flags;
    "Coefficient ring objects",          // const char *tp_doc;
    0,                                   // traverseproc tp_traverse;
    0,                                   // inquiry tp_clear;
    CoefficientRing_richcompare,         // richcmpfunc tp_richcompare;
    0,                                   // Py_ssize_t tp_weaklistoffset;
    0,                                   // getiterfunc tp_iter;
    0,                                   // iternextfunc tp_iternext;
    CoefficientRing_methods,             // struct PyMethodDef *tp_methods;
    0,                                   // struct PyMemberDef *tp_members;
    0,                                   // struct PyGetSetDef *tp_getset;
    0,                                   // struct _typeobject *tp_base;
    0,                                   // PyObject *tp_dict;
    0,                                   // descrgetfunc tp_descr_get;
    0,                                   // descrsetfunc tp_descr_set;
    0,                                   // Py_ssize_t tp_dictoffset;
    (initproc)CoefficientRing_init,      // initproc tp_init;
    0,                                   // allocfunc tp_alloc;
    CoefficientRing_new,                 // newfunc tp_new;
    0,                                   // freefunc tp_free;
    0,                                   // inquiry tp_is_gc;
    0,                                   // PyObject *tp_bases;
    0,                                   // PyObject *tp_mro;
    0,                                   // PyObject *tp_cache;
    0,                                   // PyObject *tp_subclasses;
    0,                                   // PyObject *tp_weaklist;
    0,                                   // destructor tp_del;
    0,                                   // unsigned int tp_version_tag;
    0,                                   // destructor tp_finalize;
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
	    const char* M_cstr = pythonObject2CharStar(modulus);
            lp_integer_t M;
            lp_integer_construct_from_string(lp_Z, &M, M_cstr, 10);
            int is_prime = lp_integer_is_prime(&M);
            self->K = lp_int_ring_create(&M, is_prime);
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

static PyObject *
CoefficientRing_richcompare(PyObject *self, PyObject *other, int op){
  PyObject *result = Py_NotImplemented;

  if(!PyCoefficientRing_CHECK(other)){
    //IAM: Unless I am violating some sacrosanct python contract, this seems like a no-brainer.
    if(op == Py_EQ){ return Py_False; }
    if(op == Py_NE){ return Py_True; }
  } else {
    // Get arguments
    CoefficientRing* K1 = (CoefficientRing*) self;
    CoefficientRing* K2 = (CoefficientRing*) other;

    //Are they equal?
    if (K1->K == K2->K) {
      if( op == Py_LE || op == Py_EQ || op == Py_GE){
	return Py_True;
      }
      return Py_False;
    } 

    // Is one of them Z
    if (K1->K == lp_Z) {
      // self > other
      if( op == Py_LT || op == Py_LE || op == Py_EQ){
	return Py_False;
      }
      return Py_True;
    }
    if (K2->K == lp_Z) {
      //self < other
      if( op == Py_LT || op == Py_LE){
	return Py_True;
      }
      return Py_False;
    }
    
    int cmp = lp_integer_cmp(lp_Z, &K1->K->M, &K2->K->M);

    switch (op) {
    case Py_LT:
      result = cmp < 0 ? Py_True : Py_False;
      break;
    case Py_LE:
      result = cmp <= 0 ? Py_True : Py_False;
      break;
    case Py_EQ:
      result = cmp == 0 ? Py_True : Py_False;
      break;
    case Py_NE:
      result = cmp != 0 ? Py_True : Py_False;
      break;
    case Py_GT:
      result = cmp > 0 ? Py_True : Py_False;
      break;
    case Py_GE:
      result = cmp >= 0 ? Py_True : Py_False;
      break;
    }
  }
  return result;
}
