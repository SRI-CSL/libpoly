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

#include "polypyAlgebraicNumber.h"
#include "polypyUPolynomial.h"
#include "polypyPolynomial.h"

#include <structmember.h>

static void
AlgebraicNumber_dealloc(AlgebraicNumber* self);

static PyObject*
AlgebraicNumber_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
AlgebraicNumber_init(AlgebraicNumber* self, PyObject* args);

static PyObject*
AlgebraicNumber_to_double(PyObject* self);

static PyObject*
AlgebraicNumber_richcompare(PyObject* self, PyObject* other, int op);

static PyObject*
AlgebraicNumber_str(PyObject* self);

static PyObject*
AlgebraicNumber_refine(PyObject* self);

static PyObject*
AlgebraicNumber_add(PyObject* self, PyObject* args);

static PyObject*
AlgebraicNumber_neg(PyObject* self);

static PyObject*
AlgebraicNumber_sub(PyObject* self, PyObject* args);

static PyObject*
AlgebraicNumber_mul(PyObject* self, PyObject* args);

static PyObject*
AlgebraicNumber_pow(PyObject* self, PyObject* args);

PyMethodDef AlgebraicNumber_methods[] = {
    {"refine", (PyCFunction)AlgebraicNumber_refine, METH_NOARGS, "Refines the number to half the interval"},
    {"to_double", (PyCFunction)AlgebraicNumber_to_double, METH_NOARGS, "Returns the approximation of the algebraic number"},
    {NULL}  /* Sentinel */
};

PyNumberMethods AlgebraicNumber_NumberMethods = {
     AlgebraicNumber_add, // binaryfunc nb_add;
     AlgebraicNumber_sub, // binaryfunc nb_subtract;
     AlgebraicNumber_mul, // binaryfunc nb_multiply;
     0, // binaryfunc nb_remainder;
     0, // binaryfunc nb_divmod;
     (ternaryfunc)AlgebraicNumber_pow, // ternaryfunc nb_power;
     AlgebraicNumber_neg, // unaryfunc nb_negative;
     0, // unaryfunc nb_positive;
     0, // unaryfunc nb_absolute;
     0, // inquiry nb_bool;
     0, // unaryfunc nb_invert;
     0, // binaryfunc nb_lshift;
     0, // binaryfunc nb_rshift;
     0, // binaryfunc nb_and;
     0, // binaryfunc nb_xor;
     0, // binaryfunc nb_or;
     0, // unaryfunc nb_int;
     0, // void *nb_reserved;
     0, // unaryfunc nb_float;

     0, // binaryfunc nb_inplace_add;
     0, // binaryfunc nb_inplace_subtract;
     0, // binaryfunc nb_inplace_multiply;
     0, // binaryfunc nb_inplace_remainder;
     0, // ternaryfunc nb_inplace_power;
     0, // binaryfunc nb_inplace_lshift;
     0, // binaryfunc nb_inplace_rshift;
     0, // binaryfunc nb_inplace_and;
     0, // binaryfunc nb_inplace_xor;
     0, // binaryfunc nb_inplace_or;

     0, // binaryfunc nb_floor_divide;
     0, // binaryfunc nb_true_divide;
     0, // binaryfunc nb_inplace_floor_divide;
     0, // binaryfunc nb_inplace_true_divide;

     0, // unaryfunc nb_index;

     0, // binaryfunc nb_matrix_multiply;
     0, // binaryfunc nb_inplace_matrix_multiply;
};

PyTypeObject AlgebraicNumberType = {
    PyObject_HEAD_INIT(NULL) // PyObject_VAR_HEAD
    "polypy.AlgebraicNumber", // const char *tp_name; /* For printing, in format "<module>.<name>" */
    sizeof(AlgebraicNumber),  // Py_ssize_t tp_basicsize; /* For allocation */
    0, // Py_ssize_t tp_itemsize;  /* For allocation */

    /* Methods to implement standard operations */
    (destructor)AlgebraicNumber_dealloc, //destructor tp_dealloc;
    0, //printfunc tp_print;
    0, //getattrfunc tp_getattr;
    0, //setattrfunc tp_setattr;
    0, //PyAsyncMethods *tp_as_async
    AlgebraicNumber_str, // reprfunc tp_repr;

    /* Method suites for standard classes */

    &AlgebraicNumber_NumberMethods, // PyNumberMethods *tp_as_number;
    0, //PySequenceMethods *tp_as_sequence;
    0, //PyMappingMethods *tp_as_mapping;

    /* More standard operations (here for binary compatibility) */

    0, //hashfunc tp_hash;
    0, //ternaryfunc tp_call;
    AlgebraicNumber_str, //reprfunc tp_str;
    0, //getattrofunc tp_getattro;
    0, //setattrofunc tp_setattro;

    /* Functions to access object as input/output buffer */
    0, //PyBufferProcs *tp_as_buffer;

    /* Flags to define presence of optional/expanded features */
    Py_TPFLAGS_DEFAULT, //unsigned long tp_flags;

    "Algebraic number objects", //const char *tp_doc; /* Documentation string */

    /* call function for all accessible objects */
    0, //traverseproc tp_traverse;

    /* delete references to contained objects */
    0, //inquiry tp_clear;

    /* rich comparisons */
    AlgebraicNumber_richcompare, //richcmpfunc tp_richcompare;

    /* weak reference enabler */
    0, //Py_ssize_t tp_weaklistoffset;

    /* Iterators */
    0, //getiterfunc tp_iter;
    0, //iternextfunc tp_iternext;

    /* Attribute descriptor and subclassing stuff */
    AlgebraicNumber_methods, //struct PyMethodDef *tp_methods;
    0, //struct PyMemberDef *tp_members;
    0, //struct PyGetSetDef *tp_getset;
    0, //struct _typeobject *tp_base;
    0, //PyObject *tp_dict;
    0, //descrgetfunc tp_descr_get;
    0, //descrsetfunc tp_descr_set;
    0, //Py_ssize_t tp_dictoffset;
    (initproc)AlgebraicNumber_init, //initproc tp_init;
    0, //Pallocfunc tp_alloc;
    AlgebraicNumber_new, //newfunc tp_new;
    0, //freefunc tp_free; /* Low-level free-memory routine */
    0, //inquiry tp_is_gc; /* For PyObject_IS_GC */
    0, //PyObject *tp_bases;
    0, //PyObject *tp_mro; /* method resolution order */
    0, //PyObject *tp_cache;
    0, //PyObject *tp_subclasses;
    0, //PyObject *tp_weaklist;
    0, //destructor tp_del;

    /* Type attribute cache version tag. Added in version 2.6 */
    0, //unsigned int tp_version_tag;

    0, //destructor tp_finalize;

};

static void
AlgebraicNumber_dealloc(AlgebraicNumber* self)
{
  lp_algebraic_number_destruct(&self->a);
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyAlgebraicNumber_create(const lp_algebraic_number_t* a) {
  AlgebraicNumber *self;
  self = (AlgebraicNumber*)AlgebraicNumberType.tp_alloc(&AlgebraicNumberType, 0);
  if (self != NULL) {
    if (a) {
      lp_algebraic_number_construct_copy(&self->a, a);
    } else {
      lp_dyadic_rational_t zero;
      lp_dyadic_rational_construct(&zero);
      lp_algebraic_number_construct_from_dyadic_rational(&self->a, &zero);
      lp_dyadic_rational_destruct(&zero);
    }
  }
  return (PyObject *)self;
}

static PyObject*
AlgebraicNumber_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyAlgebraicNumber_create(0);
}

/** Construct an algebraic number from a polynomial and a root index */
static int
AlgebraicNumber_init(AlgebraicNumber* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 2) {
    PyObject* f_obj = PyTuple_GetItem(args, 0);
    PyObject* root_index_obj = PyTuple_GetItem(args, 1);
    if (PyPolynomial_CHECK(f_obj) && PyLong_Check(root_index_obj)) {
      // Get the polynomial
      lp_polynomial_t* f = ((Polynomial*) f_obj)->p;
      long root_index = PyLong_AsLong(root_index_obj);
      // Get the univariate polynomial
      lp_upolynomial_t* f_u = lp_polynomial_to_univariate(f);
      if (f_u == 0) {
        // Not univariate
        return -1;
      }
      // Check the roots
      size_t roots_count = lp_upolynomial_roots_count(f_u, 0);
      if (root_index < 0 || root_index >= roots_count) {
        // Not enough roots
        lp_upolynomial_delete(f_u);
        return -1;
      }
      lp_algebraic_number_t* roots = malloc(roots_count * sizeof(lp_algebraic_number_t));
      lp_upolynomial_roots_isolate(f_u, roots, &roots_count);
      lp_algebraic_number_destruct(&self->a);
      lp_algebraic_number_construct_copy(&self->a, roots + root_index);
      int i;
      for (i = 0; i < roots_count; ++ i) {
        lp_algebraic_number_destruct(roots + i);
      }
      lp_upolynomial_delete(f_u);
      free(roots);
    } else if (PyUPolynomial_CHECK(f_obj) && PyLong_Check(root_index_obj)) {
      lp_upolynomial_t* f = ((UPolynomialObject*) f_obj)->p;
      long root_index = PyLong_AsLong(root_index_obj);
      size_t roots_count = lp_upolynomial_roots_count(f, 0);
      if (root_index < 0 || root_index >= roots_count) {
        // Not enough roots
        return -1;
      }
      lp_algebraic_number_t* roots = malloc(roots_count * sizeof(lp_algebraic_number_t));
      lp_upolynomial_roots_isolate(f, roots, &roots_count);
      lp_algebraic_number_destruct(&self->a);
      lp_algebraic_number_construct_copy(&self->a, roots + root_index);
      int i;
      for (i = 0; i < roots_count; ++i) {
        lp_algebraic_number_destruct(roots + i);
      }
      free(roots);
    } else{
      return -1;
    }
  } else {
    return -1;
  }

  // All fine, initialized
  return 0;
}

static PyObject*
AlgebraicNumber_to_double(PyObject* self) {
  AlgebraicNumber* a = (AlgebraicNumber*) self;
  return algebraic_number_to_PyFloat(&a->a);
}

static PyObject*
AlgebraicNumber_refine(PyObject* self) {
  AlgebraicNumber* a = (AlgebraicNumber*) self;
  lp_algebraic_number_refine(&a->a);
  Py_RETURN_NONE;
}

static PyObject*
AlgebraicNumber_richcompare(PyObject* self, PyObject* other, int op) {
  PyObject *result = 0;

  if (!PyAlgebraicNumber_CHECK(other)) {
    result = Py_NotImplemented;
  } else {
    lp_algebraic_number_t* self_a = &((AlgebraicNumber*) self)->a;
    lp_algebraic_number_t* other_a = &((AlgebraicNumber*) other)->a;
    int cmp = lp_algebraic_number_cmp(self_a, other_a);

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

  Py_INCREF(result);
  return result;
}

static PyObject* AlgebraicNumber_str(PyObject* self) {
  AlgebraicNumber* a = (AlgebraicNumber*) self;
  char* cstr = lp_algebraic_number_to_string(&a->a);
  PyObject* pystr = PyUnicode_FromString(cstr);
  free(cstr);
  return pystr;
}

static PyObject*
AlgebraicNumber_add(PyObject* self, PyObject* other) {

  if (PyLong_Check(self)) {
    // To enable summation
    long x = PyLong_AsLong(self);
    if (x == 0) {
      AlgebraicNumber* a = (AlgebraicNumber*) other;
      return PyAlgebraicNumber_create(&a->a);
    }
  }

  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;

  lp_algebraic_number_t sum;
  lp_algebraic_number_construct_zero(&sum);
  lp_algebraic_number_add(&sum, &a1->a, &a2->a);
  PyObject* result = PyAlgebraicNumber_create(&sum);
  lp_algebraic_number_destruct(&sum);

  return result;
}

static PyObject*
AlgebraicNumber_neg(PyObject* self) {
  if (!PyAlgebraicNumber_CHECK(self)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;

  lp_algebraic_number_t neg;
  lp_algebraic_number_construct_zero(&neg);
  lp_algebraic_number_neg(&neg, &a1->a);
  PyObject* result = PyAlgebraicNumber_create(&neg);
  lp_algebraic_number_destruct(&neg);

  return result;
}

static PyObject*
AlgebraicNumber_sub(PyObject* self, PyObject* other) {
  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;

  lp_algebraic_number_t sub;
  lp_algebraic_number_construct_zero(&sub);
  lp_algebraic_number_sub(&sub, &a1->a, &a2->a);
  PyObject* result = PyAlgebraicNumber_create(&sub);
  lp_algebraic_number_destruct(&sub);

  return result;
}

static PyObject*
AlgebraicNumber_mul(PyObject* self, PyObject* other) {
  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;

  lp_algebraic_number_t mul;
  lp_algebraic_number_construct_zero(&mul);
  lp_algebraic_number_mul(&mul, &a1->a, &a2->a);
  PyObject* result = PyAlgebraicNumber_create(&mul);
  lp_algebraic_number_destruct(&mul);

  return result;
}

static PyObject*
AlgebraicNumber_pow(PyObject* self, PyObject* other) {
  if (!PyAlgebraicNumber_CHECK(self) || !PyLong_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  long n = PyLong_AsLong(other);

  lp_algebraic_number_t pow;
  lp_algebraic_number_construct_zero(&pow);
  lp_algebraic_number_pow(&pow, &a1->a, n);
  PyObject* result = PyAlgebraicNumber_create(&pow);
  lp_algebraic_number_destruct(&pow);

  return result;
}
