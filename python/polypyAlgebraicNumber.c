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
AlgebraicNumber_add(PyObject* self, PyObject* other);

static PyObject*
AlgebraicNumber_neg(PyObject* self);

static PyObject*
AlgebraicNumber_sub(PyObject* self, PyObject* other);

static PyObject*
AlgebraicNumber_mul(PyObject* self, PyObject* other);

static PyObject*
AlgebraicNumber_div(PyObject* self, PyObject* other);

static PyObject*
AlgebraicNumber_pow(PyObject* self, PyObject* other);

static PyObject*
AlgebraicNumber_positive_root(PyObject* self, PyObject* args);

PyMethodDef AlgebraicNumber_methods[] = {
    {"refine", (PyCFunction)AlgebraicNumber_refine, METH_NOARGS, "Refines the number to half the interval"},
    {"to_double", (PyCFunction)AlgebraicNumber_to_double, METH_NOARGS, "Returns the approximation of the algebraic number"},
    {"positive_root", (PyCFunction)AlgebraicNumber_positive_root, METH_VARARGS, "Returns the positive root of the number is positive"},
    {NULL}  /* Sentinel */
};

PyNumberMethods AlgebraicNumber_NumberMethods = {
     AlgebraicNumber_add,              // binaryfunc nb_add;
     AlgebraicNumber_sub,              // binaryfunc nb_subtract;
     AlgebraicNumber_mul,              // binaryfunc nb_multiply;
     NULL,                             // binaryfunc nb_remainder;
     NULL,                             // binaryfunc nb_divmod;
     (ternaryfunc)AlgebraicNumber_pow, // ternaryfunc nb_power;
     AlgebraicNumber_neg,              // unaryfunc nb_negative;
     NULL,                             // unaryfunc nb_positive;
     NULL,                             // unaryfunc nb_absolute;
     NULL,                             // inquiry nb_bool;
     NULL,                             // unaryfunc nb_invert;
     NULL,                             // binaryfunc nb_lshift;
     NULL,                             // binaryfunc nb_rshift;
     NULL,                             // binaryfunc nb_and;
     NULL,                             // binaryfunc nb_xor;
     NULL,                             // binaryfunc nb_or;
     NULL,                             // unaryfunc nb_int;
     NULL,                             // void *nb_reserved;
     NULL,                             // unaryfunc nb_float;
     NULL,                             // binaryfunc nb_inplace_add;
     NULL,                             // binaryfunc nb_inplace_subtract;
     NULL,                             // binaryfunc nb_inplace_multiply;
     NULL,                             // binaryfunc nb_inplace_remainder;
     NULL,                             // ternaryfunc nb_inplace_power;
     NULL,                             // binaryfunc nb_inplace_lshift;
     NULL,                             // binaryfunc nb_inplace_rshift;
     NULL,                             // binaryfunc nb_inplace_and;
     NULL,                             // binaryfunc nb_inplace_xor;
     NULL,                             // binaryfunc nb_inplace_or;
     NULL,                             // binaryfunc nb_floor_divide;
     AlgebraicNumber_div,              // binaryfunc nb_true_divide;
     NULL,                             // binaryfunc nb_inplace_floor_divide;
     NULL,                             // binaryfunc nb_inplace_true_divide;
     NULL,                             // unaryfunc nb_index;
     NULL,                             // binaryfunc nb_matrix_multiply;
     NULL,                             // binaryfunc nb_inplace_matrix_multiply;
};

PyTypeObject AlgebraicNumberType = {
    {PyObject_HEAD_INIT(NULL)},          // PyObject_VAR_HEAD
    "polypy.AlgebraicNumber",            // const char *tp_name;
    sizeof(AlgebraicNumber),             // Py_ssize_t tp_basicsize; 
    0,                                   // Py_ssize_t tp_itemsize; 
    (destructor)AlgebraicNumber_dealloc, // destructor tp_dealloc;
    0,                                   // printfunc tp_print;
    NULL,                                // getattrfunc tp_getattr;
    NULL,                                // setattrfunc tp_setattr;
    NULL,                                // PyAsyncMethods *tp_as_async
    AlgebraicNumber_str,                 // reprfunc tp_repr;
    &AlgebraicNumber_NumberMethods,      // PyNumberMethods *tp_as_number;
    NULL,                                // PySequenceMethods *tp_as_sequence;
    NULL,                                // PyMappingMethods *tp_as_mapping;
    NULL,                                // hashfunc tp_hash;
    NULL,                                // ternaryfunc tp_call;
    AlgebraicNumber_str,                 // reprfunc tp_str;
    NULL,                                // getattrofunc tp_getattro;
    NULL,                                // setattrofunc tp_setattro;
    NULL,                                // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,                  // unsigned long tp_flags;
    "Algebraic number objects",          // const char *tp_doc; 
    NULL,                                // traverseproc tp_traverse;
    NULL,                                // inquiry tp_clear;
    AlgebraicNumber_richcompare,         // richcmpfunc tp_richcompare;
    0,                                   // Py_ssize_t tp_weaklistoffset;
    NULL,                                // getiterfunc tp_iter;
    NULL,                                // iternextfunc tp_iternext;
    AlgebraicNumber_methods,             // struct PyMethodDef *tp_methods;
    NULL,                                // struct PyMemberDef *tp_members;
    NULL,                                // struct PyGetSetDef *tp_getset;
    NULL,                                // struct _typeobject *tp_base;
    NULL,                                // PyObject *tp_dict;
    NULL,                                // descrgetfunc tp_descr_get;
    NULL,                                // descrsetfunc tp_descr_set;
    0,                                   // Py_ssize_t tp_dictoffset;
    (initproc)AlgebraicNumber_init,      // initproc tp_init;
    NULL,                                // Pallocfunc tp_alloc;
    AlgebraicNumber_new,                 // newfunc tp_new;
    NULL,                                // freefunc tp_free;
    NULL,                                // inquiry tp_is_gc;
    NULL,                                // PyObject *tp_bases;
    NULL,                                // PyObject *tp_mro;
    NULL,                                // PyObject *tp_cache;
    NULL,                                // PyObject *tp_subclasses;
    NULL,                                // PyObject *tp_weaklist;
    NULL,                                // destructor tp_del;
    0,                                   // unsigned int tp_version_tag;
    NULL,                                // destructor tp_finalize;
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
  return PyAlgebraicNumber_create(NULL);
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
      if (f_u == NULL) {
        // Not univariate
        return -1;
      }
      // Check the roots
      size_t roots_count = lp_upolynomial_roots_count(f_u, NULL);
      if (root_index < 0 || root_index >= (long)roots_count) {
        // Not enough roots
        lp_upolynomial_delete(f_u);
        return -1;
      }
      lp_algebraic_number_t* roots = malloc(roots_count * sizeof(lp_algebraic_number_t));
      lp_upolynomial_roots_isolate(f_u, roots, &roots_count);
      lp_algebraic_number_destruct(&self->a);
      lp_algebraic_number_construct_copy(&self->a, roots + root_index);
      for (size_t i = 0; i < roots_count; ++ i) {
        lp_algebraic_number_destruct(roots + i);
      }
      lp_upolynomial_delete(f_u);
      free(roots);
    } else if (PyUPolynomial_CHECK(f_obj) && PyLong_Check(root_index_obj)) {
      lp_upolynomial_t* f = ((UPolynomialObject*) f_obj)->p;
      long root_index = PyLong_AsLong(root_index_obj);
      size_t roots_count = lp_upolynomial_roots_count(f, NULL);
      if (root_index < 0 || root_index >= (long)roots_count) {
        // Not enough roots
        return -1;
      }
      lp_algebraic_number_t* roots = malloc(roots_count * sizeof(lp_algebraic_number_t));
      lp_upolynomial_roots_isolate(f, roots, &roots_count);
      lp_algebraic_number_destruct(&self->a);
      lp_algebraic_number_construct_copy(&self->a, roots + root_index);
      for (size_t i = 0; i < roots_count; ++i) {
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
  PyObject *result = NULL;

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
    default:
      assert(0);
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
AlgebraicNumber_div(PyObject* self, PyObject* other) {
  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;

  lp_algebraic_number_t mul;
  lp_algebraic_number_construct_zero(&mul);
  lp_algebraic_number_div(&mul, &a1->a, &a2->a);
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

static PyObject*
AlgebraicNumber_positive_root(PyObject* self, PyObject* args) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    // Get the argument
    if (PyTuple_Size(args) == 1) {
      // Get n
      PyObject* other = PyTuple_GetItem(args, 0);
      AlgebraicNumber* a1 = (AlgebraicNumber*) self;
      long n = PyLong_AsLong(other);

      lp_algebraic_number_t root;
      lp_algebraic_number_construct_zero(&root);
      lp_algebraic_number_positive_root(&root, &a1->a, n);
      PyObject* result = PyAlgebraicNumber_create(&root);
      lp_algebraic_number_destruct(&root);

      return result;
    } else {
      Py_RETURN_NONE;
    }
  } else {
    Py_RETURN_NONE;
  }
}
