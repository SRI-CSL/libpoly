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

#include "polypyUPolynomial.h"
#include "polypyInteger.h"
#include "polypyAlgebraicNumber.h"

#include "utils.h"

#include "rational_interval.h"

#include <structmember.h>

static void
UPolynomial_dealloc(UPolynomialObject* self);

static PyObject*
UPolynomial_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
UPolynomial_init(UPolynomialObject* self, PyObject* args);

static PyObject*
UPolynomial_richcompare(PyObject* self, PyObject* other, int op);

static PyObject*
UPolynomial_degree(PyObject* self);

static PyObject*
UPolynomial_coefficients(PyObject* self);

static PyObject*
UPolynomial_ring(PyObject* self);

static PyObject*
UPolynomial_gcd(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_derivative(PyObject* self);

static PyObject*
UPolynomial_extended_gcd(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_factor(PyObject* self);

static PyObject*
UPolynomial_factor_square_free(PyObject* self);

static PyObject*
UPolynomial_roots_count(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_roots_isolate(PyObject* self);

static PyObject*
UPolynomial_sturm_sequence(PyObject* self);

static PyObject*
UPolynomial_roots_find_Zp(PyObject* self);

static PyObject*
UPolynomial_evaluate(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_to_ring(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_str(PyObject* self);

static int
UPolynomialObject_nonzero(PyObject* self);

static PyObject*
UPolynomialObject_add(PyObject* self, PyObject* other);

static PyObject*
UPolynomialObject_sub(PyObject* self, PyObject* other);

static PyObject*
UPolynomialObject_mul(PyObject* self, PyObject* other);

static PyObject*
UPolynomialObject_div(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_rem(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_divmod(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_pow(PyObject* self, PyObject* other);

static PyObject*
UPolynomialObject_neg(PyObject* self);

PyMethodDef UPolynomial_methods[] = {
    {"degree", (PyCFunction)UPolynomial_degree, METH_NOARGS, "Returns the degree of the polynomial"},
    {"coefficients", (PyCFunction)UPolynomial_coefficients, METH_NOARGS, "Returns a dictionary from degrees to coefficients"},
    {"ring", (PyCFunction)UPolynomial_ring, METH_NOARGS, "Returns the ring of the polynomial"},
    {"to_ring", (PyCFunction)UPolynomial_to_ring, METH_VARARGS, "Returns the polynomial in the given ring"},
    {"gcd", (PyCFunction)UPolynomial_gcd, METH_VARARGS, "Returns the gcd of current and given polynomial in the given ring"},
    {"extended_gcd", (PyCFunction)UPolynomial_extended_gcd, METH_VARARGS, "Returns the extended gcd, i.e. (gcd, u, v), of current and given polynomial in the given ring"},
    {"factor", (PyCFunction)UPolynomial_factor, METH_NOARGS, "Returns the factorization of the polynomial"},
    {"factor_square_free", (PyCFunction)UPolynomial_factor_square_free, METH_NOARGS, "Returns the square-free factorization of the polynomial"},
    {"roots_count", (PyCFunction)UPolynomial_roots_count, METH_VARARGS, "Returns the number of real roots in the given interval"},
    {"roots_isolate", (PyCFunction)UPolynomial_roots_isolate, METH_NOARGS, "Returns the list of real roots"},
    {"sturm_sequence", (PyCFunction)UPolynomial_sturm_sequence, METH_NOARGS, "Returns the Sturm sequence"},
    {"roots_find_Zp", (PyCFunction)UPolynomial_roots_find_Zp, METH_NOARGS, "Returns the roots of the polynomial in Zp"},
    {"derivative", (PyCFunction)UPolynomial_derivative, METH_NOARGS, "Returns the derivative of the polynomial"},
    {"evaluate", (PyCFunction)UPolynomial_evaluate, METH_VARARGS, "Returns the value of the polynomial at the given point"},
    {NULL}  /* Sentinel */
};

PyNumberMethods UPolynomial_NumberMethods = {
     UPolynomialObject_add,              // binaryfunc nb_add;
     UPolynomialObject_sub,              // binaryfunc nb_subtract;
     UPolynomialObject_mul,              // binaryfunc nb_multiply;
     UPolynomialObject_rem,              // binaryfunc nb_remainder;
     UPolynomialObject_divmod,           // binaryfunc nb_divmod;
     (ternaryfunc)UPolynomialObject_pow, // ternaryfunc nb_power;
     UPolynomialObject_neg,              // unaryfunc nb_negative;
     NULL,                               // unaryfunc nb_positive;
     NULL,                               // unaryfunc nb_absolute;
     UPolynomialObject_nonzero,          // inquiry nb_bool;
     NULL,                               // unaryfunc nb_invert;
     NULL,                               // binaryfunc nb_lshift;
     NULL,                               // binaryfunc nb_rshift;
     NULL,                               // binaryfunc nb_and;
     NULL,                               // binaryfunc nb_xor;
     NULL,                               // binaryfunc nb_or;
     NULL,                               // unaryfunc nb_int;
     NULL,                               // void *nb_reserved;
     NULL,                               // unaryfunc nb_float;
     NULL,                               // binaryfunc nb_inplace_add;
     NULL,                               // binaryfunc nb_inplace_subtract;
     NULL,                               // binaryfunc nb_inplace_multiply;
     NULL,                               // binaryfunc nb_inplace_remainder;
     NULL,                               // ternaryfunc nb_inplace_power;
     NULL,                               // binaryfunc nb_inplace_lshift;
     NULL,                               // binaryfunc nb_inplace_rshift;
     NULL,                               // binaryfunc nb_inplace_and;
     NULL,                               // binaryfunc nb_inplace_xor;
     NULL,                               // binaryfunc nb_inplace_or;
     NULL,                               // binaryfunc nb_floor_divide;
     UPolynomialObject_div,              // binaryfunc nb_true_divide;
     NULL,                               // binaryfunc nb_inplace_floor_divide;
     NULL,                               // binaryfunc nb_inplace_true_divide;
     NULL,                               // unaryfunc nb_index;
     NULL,                               // binaryfunc nb_matrix_multiply;
     NULL,                               // binaryfunc nb_inplace_matrix_multiply;
};

PyTypeObject UPolynomialType = {
    {PyObject_HEAD_INIT(NULL)},      // PyObject_VAR_HEAD
    "polypy.UPolynomial",            // const char *tp_name;
    sizeof(UPolynomialObject),       // Py_ssize_t tp_basicsize;
    0,                               // Py_ssize_t tp_itemsize;
    (destructor)UPolynomial_dealloc, // destructor tp_dealloc;
    0,                               // printfunc tp_print;
    NULL,                            // getattrfunc tp_getattr;
    NULL,                            // setattrfunc tp_setattr;
    NULL,                            // PyAsyncMethods *tp_as_async;
    UPolynomial_str,                 // reprfunc tp_repr;
    &UPolynomial_NumberMethods,      // PyNumberMethods *tp_as_number;
    NULL,                            // PySequenceMethods *tp_as_sequence;
    NULL,                            // PyMappingMethods *tp_as_mapping;
    NULL,                            // hashfunc tp_hash;
    NULL,                            // ternaryfunc tp_call;
    UPolynomial_str,                 // reprfunc tp_str;
    NULL,                            // getattrofunc tp_getattro;
    NULL,                            // setattrofunc tp_setattro;
    NULL,                            // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,              // unsigned long tp_flags;
    "Univariate polynomial objects", // const char *tp_doc;
    NULL,                            // traverseproc tp_traverse;
    NULL,                            // inquiry tp_clear;
    UPolynomial_richcompare,         // richcmpfunc tp_richcompare;
    0,                               // Py_ssize_t tp_weaklistoffset;
    NULL,                            // getiterfunc tp_iter;
    NULL,                            // iternextfunc tp_iternext;
    UPolynomial_methods,             // struct PyMethodDef *tp_methods;
    NULL,                            // struct PyMemberDef *tp_members;
    NULL,                            // struct PyGetSetDef *tp_getset;
    NULL,                            // struct _typeobject *tp_base;
    NULL,                            // PyObject *tp_dict;
    NULL,                            // descrgetfunc tp_descr_get;
    NULL,                            // descrsetfunc tp_descr_set;
    0,                               // Py_ssize_t tp_dictoffset;
    (initproc)UPolynomial_init,      // initproc tp_init;
    NULL,                            // allocfunc tp_alloc;
    UPolynomial_new,                 // newfunc tp_new;
    NULL,                            // freefunc tp_free;
    NULL,                            // inquiry tp_is_gc;
    NULL,                            // PyObject *tp_bases;
    NULL,                            // PyObject *tp_mro;
    NULL,                            // PyObject *tp_cache;
    NULL,                            // PyObject *tp_subclasses;
    NULL,                            // PyObject *tp_weaklist;
    NULL,                            // destructor tp_del;
    0,                               // unsigned int tp_version_tag;
    NULL,                            // destructor tp_finalize;
};

static void
UPolynomial_dealloc(UPolynomialObject* self)
{
  if (self->p) lp_upolynomial_delete(self->p);
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyUPolynomial_create(lp_upolynomial_t* p) {
  UPolynomialObject *self;
  self = (UPolynomialObject*)UPolynomialType.tp_alloc(&UPolynomialType, 0);
  if (self != NULL) {
    self->p = p;
  }
  return (PyObject *)self;
}

static PyObject*
UPolynomial_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyUPolynomial_create(NULL);
}

static int
UPolynomial_init(UPolynomialObject* self, PyObject* args)
{
    if (PyTuple_Check(args)) {
      if (PyTuple_Size(args) == 0) {
        // Just a zero polynomial
        int coeff[] = { 0 };
        self->p = lp_upolynomial_construct_from_int(lp_Z, 0, coeff);
      } else if (PyTuple_Size(args) <= 2) {
        lp_int_ring_t* K = lp_Z;

        // Get the list of coefficients
        PyObject* coefficients = PyTuple_GetItem(args, 0);

        // The first item might be the ring
        if (PyCoefficientRing_CHECK(coefficients)) {
          K = ((CoefficientRing*) coefficients)->K;
          coefficients = PyTuple_GetItem(args, 1);
        } else {
          if (PyTuple_Size(args) == 2) {
            return -1;
          }
        }

        if (PyList_Check(coefficients)) {
          Py_ssize_t size = PyList_Size(coefficients);
          if (size > 0) {
            long c_ints[size];
            for (Py_ssize_t i = 0; i < size; ++ i) {
              PyObject* c_i = PyList_GetItem(coefficients, i);
              if (!PyLong_Check(c_i)) {
                return -1;
              } else {
                c_ints[i] = PyLong_AsLong(c_i);
              }
            }
            self->p = lp_upolynomial_construct_from_long(K, size-1, c_ints);
          } else {
            return -1;
          }
        } else {
          return -1;
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
UPolynomial_richcompare(PyObject* self, PyObject* other, int op) {
  PyObject *result = NULL;

  if (!PyUPolynomial_CHECK(other) && !PyLong_Check(other)) {
    //IAM: Surely we can be braver here?
    // if op == Py_EQ we could return Py_False
    // if op == Py_NE we could return Py_True
    // I am sure this holds for other versions of:
    // _cmp AND _richcompare in this codebase.
    result = Py_NotImplemented;
  } else {
    lp_upolynomial_t* self_p = ((UPolynomialObject*) self)->p;
    lp_upolynomial_t* other_p = NULL;
    if (PyUPolynomial_CHECK(other)) {
      other_p = ((UPolynomialObject*) other)->p;
    } else {
      long c = PyLong_AsLong(other);
      const lp_int_ring_t* K = lp_upolynomial_ring(self_p);
      other_p = lp_upolynomial_construct_from_long(K, 0, &c);
    }

    int cmp = lp_upolynomial_cmp(self_p, other_p);

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

    if (PyLong_Check(other)) {
      lp_upolynomial_delete(other_p);
    }
  }

  Py_INCREF(result);
  return result;
}

static PyObject*
UPolynomial_degree(PyObject* self) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    return PyLong_FromUnsignedLong(lp_upolynomial_degree(p->p));
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject*
UPolynomial_coefficients(PyObject* self) {
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  size_t size = lp_upolynomial_degree(p) + 1;

  lp_integer_t coefficients[size];
  for (size_t i = 0; i < size; ++ i) {
    lp_integer_construct_from_int(lp_Z, coefficients + i, 0);
  }

  lp_upolynomial_unpack(p, coefficients);

  PyObject* list = PyList_New(size);

  for (size_t i = 0; i < size; ++ i) {
    PyObject* c = integer_to_PyLong(coefficients + i);
    PyList_SetItem(list, i, c);
    lp_integer_destruct(coefficients + i);
  }

  return list;
}

static PyObject*
UPolynomial_ring(PyObject* self) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    lp_int_ring_t* K = (lp_int_ring_t*) lp_upolynomial_ring(p->p);
    lp_int_ring_attach(K);
    return PyCoefficientRing_create(K);
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject*
UPolynomial_to_ring(PyObject* self, PyObject* args) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    // Get the ring argument
    if (PyTuple_Size(args) == 1) {
      // Get the ring
      PyObject* arg = PyTuple_GetItem(args, 0);
      if (PyCoefficientRing_CHECK(arg)) {
        CoefficientRing* K = (CoefficientRing*)arg;
        lp_upolynomial_t* p_K = lp_upolynomial_construct_copy_K(K->K, p->p);
        return PyUPolynomial_create(p_K);
      } else {
        Py_RETURN_NONE;
      }
    } else {
      Py_RETURN_NONE;
    }
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject* UPolynomial_str(PyObject* self) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    char* p_str = lp_upolynomial_to_string(p->p);
    PyObject* str = PyUnicode_FromString(p_str);
    free(p_str);
    return str;
  } else {
    Py_RETURN_NONE;
  }
}

// Add other integer
static PyObject*
UPolynomialObject_add_number(PyObject* self, PyObject* other) {
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  lp_integer_t c;
  const lp_int_ring_t* K = lp_upolynomial_ring(p1->p);
  PyLong_or_Int_to_integer(other, K, &c);
  lp_upolynomial_t* c_p = lp_upolynomial_construct(K, 0, &c);
  lp_upolynomial_t* sum = lp_upolynomial_add(p1->p, c_p);
  lp_upolynomial_delete(c_p);
  lp_integer_destruct(&c);
  return PyUPolynomial_create(sum);
}

static PyObject*
UPolynomialObject_add(PyObject* self, PyObject* other) {
  // Integer addition
  if (PyLong_or_Int_Check(other)) {
    return UPolynomialObject_add_number(self, other);
  }
  if (PyLong_or_Int_Check(self)) {
      return UPolynomialObject_add_number(other, self);
  }
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) other;
  // Add the polynomials
  lp_upolynomial_t* sum = lp_upolynomial_add(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(sum);
}

static PyObject*
UPolynomialObject_sub_int(PyObject* self, PyObject* other, int negate) {
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  lp_integer_t c;
  const lp_int_ring_t* K = lp_upolynomial_ring(p1->p);
  PyLong_or_Int_to_integer(other, K, &c);
  lp_upolynomial_t* c_p = lp_upolynomial_construct(K, 0, &c);
  lp_upolynomial_t* sub =
      negate ? lp_upolynomial_sub(c_p, p1->p) : lp_upolynomial_sub(p1->p, c_p);
  lp_upolynomial_delete(c_p);
  lp_integer_destruct(&c);
  return PyUPolynomial_create(sub);
}

static PyObject*
UPolynomialObject_sub(PyObject* self, PyObject* other) {
  // Integer addition
  if (PyLong_or_Int_Check(other)) {
    return UPolynomialObject_sub_int(self, other, 0);
  }
  if (PyLong_or_Int_Check(self)) {
    return UPolynomialObject_sub_int(other, self, 1);
  }
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) other;
  // Add the polynomials
  lp_upolynomial_t* sum = lp_upolynomial_sub(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(sum);
}

// Multiply other integer
static PyObject*
UPolynomialObject_mul_int(PyObject* self, PyObject* other) {
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  lp_integer_t c;
  PyLong_or_Int_to_integer(other, lp_upolynomial_ring(p1->p), &c);
  lp_upolynomial_t* sum = lp_upolynomial_mul_c(p1->p, &c);
  lp_integer_destruct(&c);
  return PyUPolynomial_create(sum);
}

static PyObject*
UPolynomialObject_mul(PyObject* self, PyObject* other) {
  // Integer multiplication
  if (PyLong_or_Int_Check(other)) {
    return UPolynomialObject_mul_int(self, other);
  }
  if (PyLong_or_Int_Check(self)) {
    return UPolynomialObject_mul_int(other, self);
  }
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) other;
  // Add the polynomials
  lp_upolynomial_t* sum = lp_upolynomial_mul(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(sum);
}

static PyObject*
UPolynomialObject_div(PyObject* self, PyObject* args) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(args)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) args;
  // Divide the polynomials
  lp_upolynomial_t* div = lp_upolynomial_div_exact(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(div);
}

static PyObject*
UPolynomialObject_rem(PyObject* self, PyObject* args) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(args)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) args;
  // Divide the polynomials
  lp_upolynomial_t* rem = lp_upolynomial_rem_exact(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(rem);
}

static PyObject*
UPolynomialObject_divmod(PyObject* self, PyObject* args) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(args)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) args;
  // To store div, rem
  lp_upolynomial_t* div = NULL;
  lp_upolynomial_t* rem = NULL;
  // Divide the polynomials
  lp_upolynomial_div_rem_exact(p1->p, p2->p, &div, &rem);
  // Return the result
  PyObject* pair = PyTuple_New(2);
  PyObject* divObj = PyUPolynomial_create(div);
  PyObject* remObj = PyUPolynomial_create(rem);
  Py_INCREF(divObj);
  Py_INCREF(remObj);
  PyTuple_SetItem(pair, 0, divObj);
  PyTuple_SetItem(pair, 1, remObj);
  return pair;
}

// object.__neg__(self)
static PyObject*
UPolynomialObject_neg(PyObject* self) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Add the polynomials
  lp_integer_t c;
  lp_integer_construct_from_int(lp_upolynomial_ring(p->p), &c, -1);
  lp_upolynomial_t* neg = lp_upolynomial_mul_c(p->p, &c);
  lp_integer_destruct(&c);
  // Return the result
  return PyUPolynomial_create(neg);
}

static PyObject*
UPolynomialObject_pow(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyLong_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  long pow = PyLong_AsLong(other);
  // Power the polynomial
  lp_upolynomial_t* p_pow = lp_upolynomial_pow(p->p, pow);
  // Return the result
  return PyUPolynomial_create(p_pow);
}

static int
UPolynomialObject_nonzero(PyObject* self) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Return the result
  return !lp_upolynomial_is_zero(p->p);
}

static PyObject*
UPolynomial_gcd(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* arg = PyTuple_GetItem(args, 0);
    if (PyUPolynomial_CHECK(arg)) {
      lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
      lp_upolynomial_t* q = ((UPolynomialObject*) arg)->p;
      lp_upolynomial_t* gcd = lp_upolynomial_gcd(p, q);
      return PyUPolynomial_create(gcd);
    } else {
      Py_RETURN_NONE;
    }
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject*
UPolynomial_extended_gcd(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* arg = PyTuple_GetItem(args, 0);
    if (PyUPolynomial_CHECK(arg)) {
      lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
      lp_upolynomial_t* q = ((UPolynomialObject*) arg)->p;
      lp_upolynomial_t* u = NULL;
      lp_upolynomial_t* v = NULL;
      lp_upolynomial_t* gcd = lp_upolynomial_extended_gcd(p, q, &u, &v);

      PyObject* t = PyTuple_New(3);
      PyObject* t0 = PyUPolynomial_create(gcd);
      PyObject* t1 = PyUPolynomial_create(u);
      PyObject* t2 = PyUPolynomial_create(v);
      Py_INCREF(t0);
      Py_INCREF(t1);
      Py_INCREF(t2);
      PyTuple_SetItem(t, 0, t0);
      PyTuple_SetItem(t, 1, t1);
      PyTuple_SetItem(t, 2, t2);
      return t;
    } else {
      Py_RETURN_NONE;
    }
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject* factors_to_PyList(lp_upolynomial_factors_t* factors) {
  // Construct the result
  size_t size = lp_upolynomial_factors_size(factors);
  PyObject* factors_list = PyList_New(size + 1);

  // Copy the constant
  PyObject* constant = integer_to_PyLong(lp_upolynomial_factors_get_constant(factors));
  PyList_SetItem(factors_list, 0, constant); // Steals the reference

  // Copy over the factors
  for (size_t i = 0; i < size; ++ i) {
    size_t degree;
    PyObject* p_i = PyUPolynomial_create(lp_upolynomial_factors_get_factor(factors, i, &degree));
    Py_INCREF(p_i);
    PyObject* d = PyLong_FromSize_t(degree);
    PyObject* pair = PyTuple_New(2);
    PyTuple_SetItem(pair, 0, p_i);
    PyTuple_SetItem(pair, 1, d);
    PyList_SetItem(factors_list, i + 1, pair);
  }

  // Return the list
  return factors_list;
}

static PyObject*
UPolynomial_factor(PyObject* self) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Factor
  lp_upolynomial_factors_t* factors = lp_upolynomial_factor(p->p);
  // Create the list
  PyObject* factors_list = factors_to_PyList(factors);
  // Get rid of the factors (not the polynomials)
  lp_upolynomial_factors_destruct(factors, 0);
  // Return the list
  return factors_list;
}

static PyObject*
UPolynomial_factor_square_free(PyObject* self) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Factor
  lp_upolynomial_factors_t* factors = lp_upolynomial_factor_square_free(p->p);
  // Create the list
  PyObject* factors_list = factors_to_PyList(factors);
  // Get rid of the factors (not the polynomials)
  lp_upolynomial_factors_destruct(factors, 0);
  // Return the list
  return factors_list;
}

static PyObject*
UPolynomial_roots_count(PyObject* self, PyObject* args) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;

  int roots = 0;

  if (PyTuple_Check(args) && PyTuple_Size(args) == 2) {
    PyObject* a = PyTuple_GetItem(args, 0);
    PyObject* b = PyTuple_GetItem(args, 1);

    // Ends of the interval
    lp_dyadic_rational_t a_rat, b_rat;

    if (PyLong_or_Int_Check(a)) {
      lp_integer_t a_int;
      PyLong_or_Int_to_integer(a, lp_Z, &a_int);
      lp_dyadic_rational_construct_from_integer(&a_rat, &a_int);
      lp_integer_destruct(&a_int);
    } else if (PyFloat_Check(a)) {
      PyFloat_to_dyadic_rational(a, &a_rat);
    } else {
      Py_RETURN_NONE;
    }

    if (PyLong_or_Int_Check(b)) {
      lp_integer_t b_int;
      PyLong_or_Int_to_integer(b, lp_Z, &b_int);
      lp_dyadic_rational_construct_from_integer(&b_rat, &b_int);
      lp_integer_destruct(&b_int);
    } else if (PyFloat_Check(b)) {
      PyFloat_to_dyadic_rational(b, &b_rat);
    } else {
      lp_dyadic_rational_destruct(&a_rat);
      Py_RETURN_NONE;
    }

    // The interval
    lp_rational_interval_t ab;
    lp_rational_interval_construct_from_dyadic(&ab, &a_rat, 1, &b_rat, 1);

    // Count
    roots = lp_upolynomial_roots_count(p->p, &ab);

    // Remove the temporaries
    lp_rational_interval_destruct(&ab);
    lp_dyadic_rational_destruct(&a_rat);
    lp_dyadic_rational_destruct(&b_rat);

  } else if (PyTuple_Size(args) == 0) {
    // count in (-inf, inf)
    roots = lp_upolynomial_roots_count(p->p, NULL);
  } else {
    Py_RETURN_NONE;
  }

  return PyLong_FromLong(roots);
}

static PyObject*
UPolynomial_roots_isolate(PyObject* self) {

  // Get the polynomial
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;

  // The isolating intervals
  size_t roots_size = lp_upolynomial_degree(p)+1;
  lp_algebraic_number_t* roots = malloc(sizeof(lp_algebraic_number_t)*roots_size);

  // Isolate the intervals (up to 2^precision)
  lp_upolynomial_roots_isolate(p, roots, &roots_size);

  // Generate a list of floats
  PyObject* list = PyList_New(roots_size);

  for (size_t i = 0; i < roots_size; ++ i) {
    PyObject* c = PyAlgebraicNumber_create(&roots[i]);
    PyList_SetItem(list, i, c);
  }

  // Get rid of the temporaries
  for (size_t i = 0; i < roots_size; ++ i) {
    lp_algebraic_number_destruct(roots + i);
  }
  free(roots);

  // This is the list
  return list;
}

static PyObject*
UPolynomial_derivative(PyObject* self) {
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  lp_upolynomial_t* p_derivative = lp_upolynomial_derivative(p);
  return PyUPolynomial_create(p_derivative);
}

static PyObject* upolynomials_to_PyList(lp_upolynomial_t** list, size_t size) {
  // Construct the result
  PyObject* pylist = PyList_New(size);

  // Copy over the polynomials
  for (size_t i = 0; i < size; ++ i) {
    PyObject* pylist_i = PyUPolynomial_create(list[i]);
    Py_INCREF(pylist_i);
    PyList_SetItem(pylist, i, pylist_i);
  }

  // Return the list
  return pylist;
}


static PyObject*
UPolynomial_sturm_sequence(PyObject* self) {
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  lp_upolynomial_t** S;
  size_t S_size;
  lp_upolynomial_sturm_sequence(p, &S, &S_size);
  PyObject* result = upolynomials_to_PyList(S, S_size);
  free(S);
  return result;
}

static PyObject* integer_list_to_PyList(const lp_integer_t *list, size_t size) {
  // Construct the result
  PyObject* pylist = PyList_New(size);

  // Copy over the integers
  for (size_t i = 0; i < size; ++ i) {
    PyObject* pylist_i = integer_to_PyLong(&list[i]);
    Py_INCREF(pylist_i);
    PyList_SetItem(pylist, i, pylist_i);
  }

  // Return the list
  return pylist;
}

static PyObject*
UPolynomial_roots_find_Zp(PyObject* self) {
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  lp_integer_t* roots;
  size_t roots_size;
  lp_upolynomial_roots_find_Zp(p, &roots, &roots_size);
  PyObject* result = integer_list_to_PyList(roots, roots_size);
  for (size_t i = 0; i < roots_size; ++i) {
    lp_integer_destruct(&roots[i]);
  }
  free(roots);
  return result;
}

static PyObject*
UPolynomial_evaluate(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
      PyObject* x = PyTuple_GetItem(args, 0);
      // Ends of the interval
      lp_dyadic_rational_t x_rat;
      if (PyLong_or_Int_Check(x)) {
        lp_integer_t x_int;
        PyLong_or_Int_to_integer(x, lp_Z, &x_int);
        lp_dyadic_rational_construct_from_integer(&x_rat, &x_int);
        lp_integer_destruct(&x_int);
      } else if (PyFloat_Check(x)) {
        PyFloat_to_dyadic_rational(x, &x_rat);
      } else {
        Py_RETURN_NONE;
      }

      lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
      lp_dyadic_rational_t value;
      lp_dyadic_rational_construct(&value);
      lp_upolynomial_evaluate_at_dyadic_rational(p, &x_rat, &value);
      PyObject* result = dyadic_rational_to_PyFloat(&value);
      lp_dyadic_rational_destruct(&x_rat);
      lp_dyadic_rational_destruct(&value);
      return result;
  } else {
    Py_RETURN_NONE;
  }

}
