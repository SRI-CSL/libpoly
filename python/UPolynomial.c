/*
 * UPolynomial.c
 *
 *  Created on: Nov 8, 2013
 *      Author: dejan
 */

#include "UPolynomial.h"
#include "Integer.h"
#include "AlgebraicNumber.h"

#include "utils.h"

#include <structmember.h>

static void
UPolynomial_dealloc(UPolynomialObject* self);

static PyObject*
UPolynomial_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
UPolynomial_init(UPolynomialObject* self, PyObject* args);

static int
UPolynomial_cmp(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_richcompare(PyObject* self, PyObject* args, int op);

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
UPolynomial_evaluate(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_to_ring(PyObject* self, PyObject* args);

static PyObject*
UPolynomial_str(PyObject* self);

static int
UPolynomialObject_nonzero(PyObject* self);

static PyObject*
UPolynomialObject_add(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_sub(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_mul(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_div(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_rem(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_divmod(PyObject* self, PyObject* args);

static PyObject*
UPolynomialObject_pow(PyObject* self, PyObject* args);

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
    {"derivative", (PyCFunction)UPolynomial_derivative, METH_NOARGS, "Returns the derivative of the polynomial"},
    {"evaluate", (PyCFunction)UPolynomial_evaluate, METH_VARARGS, "Returns the value of the polynomial at the given point"},
    {NULL}  /* Sentinel */
};

PyNumberMethods UPolynomial_NumberMethods = {
     UPolynomialObject_add, // binaryfunc nb_add;
     UPolynomialObject_sub, // binaryfunc nb_subtract;
     UPolynomialObject_mul, // binaryfunc nb_multiply;
     UPolynomialObject_div, // binaryfunc nb_divide;
     UPolynomialObject_rem, // binaryfunc nb_remainder;
     UPolynomialObject_divmod, // binaryfunc nb_divmod;
     (ternaryfunc)UPolynomialObject_pow, // ternaryfunc nb_power;
     UPolynomialObject_neg, // unaryfunc nb_negative;
     0, // unaryfunc nb_positive;
     0, // unaryfunc nb_absolute;
     UPolynomialObject_nonzero, // inquiry nb_nonzero;       /* Used by PyObject_IsTrue */
     0, // unaryfunc nb_invert;
     0, // binaryfunc nb_lshift;
     0, // binaryfunc nb_rshift;
     0, // binaryfunc nb_and;
     0, // binaryfunc nb_xor;
     0, // binaryfunc nb_or;
     0, // coercion nb_coerce;       /* Used by the coerce() function */
     0, // unaryfunc nb_int;
     0, // unaryfunc nb_long;
     0, // unaryfunc nb_float;
     0, // unaryfunc nb_oct;
     0, // unaryfunc nb_hex;

     /* Added in release 2.0 */
     0, // binaryfunc nb_inplace_add;
     0, // binaryfunc nb_inplace_subtract;
     0, // binaryfunc nb_inplace_multiply;
     0, // binaryfunc nb_inplace_divide;
     0, // binaryfunc nb_inplace_remainder;
     0, // ternaryfunc nb_inplace_power;
     0, // binaryfunc nb_inplace_lshift;
     0, // binaryfunc nb_inplace_rshift;
     0, // binaryfunc nb_inplace_and;
     0, // binaryfunc nb_inplace_xor;
     0, // binaryfunc nb_inplace_or;

     /* Added in release 2.2 */
     0, // binaryfunc nb_floor_divide;
     0, // binaryfunc nb_true_divide;
     0, // binaryfunc nb_inplace_floor_divide;
     0, // binaryfunc nb_inplace_true_divide;

     /* Added in release 2.5 */
     0 // unaryfunc nb_index;
};

PyTypeObject UPolynomialType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "polypy.UPolynomial",       /*tp_name*/
    sizeof(UPolynomialObject),  /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)UPolynomial_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    UPolynomial_cmp,      /*tp_compare*/
    UPolynomial_str,            /*tp_repr*/
    &UPolynomial_NumberMethods, /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    UPolynomial_str,            /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "Univariate polynomial objects", /* tp_doc */
    0,                             /* tp_traverse */
    0,                         /* tp_clear */
    UPolynomial_richcompare,   /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    UPolynomial_methods,       /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)UPolynomial_init,/* tp_init */
    0,                         /* tp_alloc */
    UPolynomial_new,           /* tp_new */
};

static void
UPolynomial_dealloc(UPolynomialObject* self)
{
  if (self->p) upolynomial_ops.delete(self->p);
  self->ob_type->tp_free((PyObject*)self);
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
  return PyUPolynomial_create(0);
}

static int
UPolynomial_init(UPolynomialObject* self, PyObject* args)
{
    if (PyTuple_Check(args)) {
      if (PyTuple_Size(args) == 0) {
        // Just a zero polynomial
        int coeff[] = { 0 };
        self->p = upolynomial_ops.construct_from_int(lp_Z, 0, coeff);
      } else if (PyTuple_Size(args) <= 2) {
        lp_int_ring K = lp_Z;

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
            Py_ssize_t i;
            for (i = 0; i < size; ++ i) {
              PyObject* c_i = PyList_GetItem(coefficients, i);
              if (!PyInt_Check(c_i)) {
                return -1;
              } else {
                c_ints[i] = PyInt_AsLong(c_i);
              }
            }
            self->p = upolynomial_ops.construct_from_long(K, size-1, c_ints);
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
  PyObject *result = 0;

  if (!PyUPolynomial_CHECK(other) || !PyInt_Check(other)) {
    result = Py_NotImplemented;
  } else {
    lp_upolynomial_t* self_p = ((UPolynomialObject*) self)->p;
    lp_upolynomial_t* other_p = 0;
    if (PyUPolynomial_CHECK(other)) {
      other_p = ((UPolynomialObject*) other)->p;
    } else {
      long c = PyInt_AsLong(other);
      lp_int_ring K = upolynomial_ops.ring(self_p);
      other_p = upolynomial_ops.construct_from_long(K, 0, &c);
    }

    int cmp = upolynomial_ops.cmp(self_p, other_p);

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

    if (PyInt_Check(other)) {
      upolynomial_ops.delete(other_p);
    }
  }

  Py_INCREF(result);
  return result;
}

static int
UPolynomial_cmp(PyObject* self, PyObject* other) {

  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyUPolynomial_CHECK(other)) {
    // should return -1 and set an exception condition when an error occurred
    return -1;
  }
  // Get arguments
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  UPolynomialObject* p2 = (UPolynomialObject*) other;
  // Compare
  int cmp = upolynomial_ops.cmp(p1->p, p2->p);
  return cmp > 0 ? 1 : cmp < 0 ? -1 : 0;
}

static PyObject*
UPolynomial_degree(PyObject* self) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    return PyInt_FromLong(upolynomial_ops.degree(p->p));
  } else {
    Py_RETURN_NONE;
  }
}

static PyObject*
UPolynomial_coefficients(PyObject* self) {
  int i;

  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  size_t size = upolynomial_ops.degree(p) + 1;

  lp_integer_t coefficients[size];
  for (i = 0; i < size; ++ i) {
    lp_integer_ops.construct_from_int(lp_Z, coefficients + i, 0);
  }

  upolynomial_ops.unpack(p, coefficients);

  PyObject* list = PyList_New(size);

  for (i = 0; i < size; ++ i) {
    PyObject* c = integer_to_PyInt(coefficients + i);
    PyList_SetItem(list, i, c);
    lp_integer_ops.destruct(coefficients + i);
  }

  return list;
}

static PyObject*
UPolynomial_ring(PyObject* self) {
  UPolynomialObject* p = (UPolynomialObject*) self;
  if (p) {
    lp_int_ring K = upolynomial_ops.ring(p->p);
    lp_int_ring_ops.attach(K);
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
        lp_upolynomial_t* p_K = upolynomial_ops.construct_copy_K(K->K, p->p);
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
    char* p_str = upolynomial_ops.to_string(p->p);
    PyObject* str = PyString_FromString(p_str);
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
  lp_int_ring K = upolynomial_ops.ring(p1->p);
  PyLong_or_Int_to_integer(other, K, &c);
  lp_upolynomial_t* c_p = upolynomial_ops.construct(K, 0, &c);
  lp_upolynomial_t* sum = upolynomial_ops.add(p1->p, c_p);
  upolynomial_ops.delete(c_p);
  lp_integer_ops.destruct(&c);
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
  lp_upolynomial_t* sum = upolynomial_ops.add(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(sum);
}

static PyObject*
UPolynomialObject_sub_int(PyObject* self, PyObject* other, int negate) {
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  lp_integer_t c;
  lp_int_ring K = upolynomial_ops.ring(p1->p);
  PyLong_or_Int_to_integer(other, K, &c);
  lp_upolynomial_t* c_p = upolynomial_ops.construct(K, 0, &c);
  lp_upolynomial_t* sub =
      negate ? upolynomial_ops.sub(c_p, p1->p) : upolynomial_ops.sub(p1->p, c_p);
  upolynomial_ops.delete(c_p);
  lp_integer_ops.destruct(&c);
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
  lp_upolynomial_t* sum = upolynomial_ops.sub(p1->p, p2->p);
  // Return the result
  return PyUPolynomial_create(sum);
}

// Multiply other integer
static PyObject*
UPolynomialObject_mul_int(PyObject* self, PyObject* other) {
  UPolynomialObject* p1 = (UPolynomialObject*) self;
  lp_integer_t c;
  PyLong_or_Int_to_integer(other, upolynomial_ops.ring(p1->p), &c);
  lp_upolynomial_t* sum = upolynomial_ops.multiply_c(p1->p, &c);
  lp_integer_ops.destruct(&c);
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
  lp_upolynomial_t* sum = upolynomial_ops.multiply(p1->p, p2->p);
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
  lp_upolynomial_t* div = upolynomial_ops.div_exact(p1->p, p2->p);
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
  lp_upolynomial_t* rem = upolynomial_ops.rem_exact(p1->p, p2->p);
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
  lp_upolynomial_t* div = 0;
  lp_upolynomial_t* rem = 0;
  // Divide the polynomials
  upolynomial_ops.div_rem_exact(p1->p, p2->p, &div, &rem);
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
  lp_integer_ops.construct_from_int(upolynomial_ops.ring(p->p), &c, -1);
  lp_upolynomial_t* neg = upolynomial_ops.multiply_c(p->p, &c);
  lp_integer_ops.destruct(&c);
  // Return the result
  return PyUPolynomial_create(neg);
}

static PyObject*
UPolynomialObject_pow(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyUPolynomial_CHECK(self) || !PyInt_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  long pow = PyInt_AsLong(other);
  // Power the polynomial
  lp_upolynomial_t* p_pow = upolynomial_ops.power(p->p, pow);
  // Return the result
  return PyUPolynomial_create(p_pow);
}

static int
UPolynomialObject_nonzero(PyObject* self) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Return the result
  return !upolynomial_ops.is_zero(p->p);
}

static PyObject*
UPolynomial_gcd(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* arg = PyTuple_GetItem(args, 0);
    if (PyUPolynomial_CHECK(arg)) {
      lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
      lp_upolynomial_t* q = ((UPolynomialObject*) arg)->p;
      lp_upolynomial_t* gcd = upolynomial_ops.gcd(p, q);
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
      lp_upolynomial_t* u = 0;
      lp_upolynomial_t* v = 0;
      lp_upolynomial_t* gcd = upolynomial_ops.extended_gcd(p, q, &u, &v);

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
  size_t size = lp_upolynomial_factors_ops.size(factors);
  PyObject* factors_list = PyList_New(size + 1);

  // Copy the constant
  PyObject* constant = integer_to_PyInt(lp_upolynomial_factors_ops.get_constant(factors));
  PyList_SetItem(factors_list, 0, constant); // Steals the reference

  // Copy over the factors
  int i;
  for (i = 0; i < size; ++ i) {
    size_t degree;
    PyObject* p_i = PyUPolynomial_create(lp_upolynomial_factors_ops.get_factor(factors, i, &degree));
    Py_INCREF(p_i);
    PyObject* d = PyInt_FromSize_t(degree);
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
  lp_upolynomial_factors_t* factors = upolynomial_ops.factor(p->p);
  // Create the list
  PyObject* factors_list = factors_to_PyList(factors);
  // Get rid of the factors (not the polynomials)
  lp_upolynomial_factors_ops.destruct(factors, 0);
  // Return the list
  return factors_list;
}

static PyObject*
UPolynomial_factor_square_free(PyObject* self) {
  // Get arguments
  UPolynomialObject* p = (UPolynomialObject*) self;
  // Factor
  lp_upolynomial_factors_t* factors = upolynomial_ops.factor_square_free(p->p);
  // Create the list
  PyObject* factors_list = factors_to_PyList(factors);
  // Get rid of the factors (not the polynomials)
  lp_upolynomial_factors_ops.destruct(factors, 0);
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
      lp_dyadic_rational_ops.construct_from_integer(&a_rat, &a_int);
      lp_integer_ops.destruct(&a_int);
    } else if (PyFloat_Check(a)) {
      PyFloat_to_dyadic_rational(a, &a_rat);
    } else {
      Py_RETURN_NONE;
    }

    if (PyLong_or_Int_Check(b)) {
      lp_integer_t b_int;
      PyLong_or_Int_to_integer(b, lp_Z, &b_int);
      lp_dyadic_rational_ops.construct_from_integer(&b_rat, &b_int);
      lp_integer_ops.destruct(&b_int);
    } else if (PyFloat_Check(b)) {
      PyFloat_to_dyadic_rational(b, &b_rat);
    } else {
      lp_dyadic_rational_ops.destruct(&a_rat);
      Py_RETURN_NONE;
    }

    // The interval
    interval_t ab;
    interval_ops.construct_from_dyadic(&ab, &a_rat, 1, &b_rat, 1);

    // Count
    roots = upolynomial_ops.roots_count(p->p, &ab);

    // Remove the temporaries
    interval_ops.destruct(&ab);
    lp_dyadic_rational_ops.destruct(&a_rat);
    lp_dyadic_rational_ops.destruct(&b_rat);

  } else if (PyTuple_Size(args) == 0) {
    // count in (-inf, inf)
    roots = upolynomial_ops.roots_count(p->p, 0);
  } else {
    Py_RETURN_NONE;
  }

  return PyInt_FromLong(roots);
}

static PyObject*
UPolynomial_roots_isolate(PyObject* self) {

  // Get the polynomial
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;

  // The isolating intervals
  size_t roots_size = upolynomial_ops.degree(p)+1;
  lp_algebraic_number_t* roots = malloc(sizeof(lp_algebraic_number_t)*roots_size);

  // Isolate the intervals (up to 2^precision)
  upolynomial_ops.roots_isolate(p, roots, &roots_size);

  // Generate a list of floats
  PyObject* list = PyList_New(roots_size);

  int i;
  for (i = 0; i < roots_size; ++ i) {
    PyObject* c = PyAlgebraicNumber_create(&roots[i]);
    PyList_SetItem(list, i, c);
  }

  // Get rid of the temporaries
  for (i = 0; i < roots_size; ++ i) {
    lp_algebraic_number_destruct(roots + i);
  }
  free(roots);

  // This is the list
  return list;
}

static PyObject*
UPolynomial_derivative(PyObject* self) {
  lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
  lp_upolynomial_t* p_derivative = upolynomial_ops.derivative(p);
  return PyUPolynomial_create(p_derivative);
}

static PyObject* upolynomials_to_PyList(lp_upolynomial_t** list, size_t size) {
  // Construct the result
  PyObject* pylist = PyList_New(size);

  // Copy over the polynomials
  size_t i;
  for (i = 0; i < size; ++ i) {
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
  upolynomial_ops.sturm_sequence(p, &S, &S_size);
  PyObject* result = upolynomials_to_PyList(S, S_size);
  free(S);
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
        lp_dyadic_rational_ops.construct_from_integer(&x_rat, &x_int);
        lp_integer_ops.destruct(&x_int);
      } else if (PyFloat_Check(x)) {
        PyFloat_to_dyadic_rational(x, &x_rat);
      } else {
        Py_RETURN_NONE;
      }

      lp_upolynomial_t* p = ((UPolynomialObject*) self)->p;
      lp_dyadic_rational_t value;
      lp_dyadic_rational_ops.construct(&value);
      upolynomial_ops.evaluate_at_dyadic_rational(p, &x_rat, &value);
      PyObject* result = dyadic_rational_to_PyFloat(&value);
      lp_dyadic_rational_ops.destruct(&x_rat);
      lp_dyadic_rational_ops.destruct(&value);
      return result;
  } else {
    Py_RETURN_NONE;
  }

}
