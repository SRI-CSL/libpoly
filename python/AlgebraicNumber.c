/*
 * AlgebraicNumber.c
 *
 *  Created on: Jan 26, 2014
 *      Author: dejan
 */

#include "AlgebraicNumber.h"
#include "UPolynomial.h"

#include <structmember.h>

static void
AlgebraicNumber_dealloc(AlgebraicNumber* self);

static PyObject*
AlgebraicNumber_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
AlgebraicNumber_init(AlgebraicNumber* self, PyObject* args);

static PyObject*
AlgebraicNumber_to_double(PyObject* self);

  static int
AlgebraicNumber_cmp(PyObject* self, PyObject* args);

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
     0, // binaryfunc nb_divide;
     0, // binaryfunc nb_remainder;
     0, // binaryfunc nb_divmod;
     (ternaryfunc)AlgebraicNumber_pow, // ternaryfunc nb_power;
     AlgebraicNumber_neg, // unaryfunc nb_negative;
     0, // unaryfunc nb_positive;
     0, // unaryfunc nb_absolute;
     0, // inquiry nb_nonzero;       /* Used by PyObject_IsTrue */
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

PyTypeObject AlgebraicNumberType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "polypy.AlgebraicNumber",   /*tp_name*/
    sizeof(AlgebraicNumber),    /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)AlgebraicNumber_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    AlgebraicNumber_cmp,        /*tp_compare*/
    AlgebraicNumber_str,        /*tp_repr*/
    &AlgebraicNumber_NumberMethods, /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    AlgebraicNumber_str,        /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "Algebraic number objects", /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    AlgebraicNumber_richcompare,/* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    AlgebraicNumber_methods,    /* tp_methods */
    0,                          /* tp_members */
    0,                          /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)AlgebraicNumber_init,/* tp_init */
    0,                          /* tp_alloc */
    AlgebraicNumber_new,        /* tp_new */
};

static void
AlgebraicNumber_dealloc(AlgebraicNumber* self)
{
  algebraic_number_ops.destruct(&self->a);
  self->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyAlgebraicNumber_create(const algebraic_number_t* a) {
  AlgebraicNumber *self;
  self = (AlgebraicNumber*)AlgebraicNumberType.tp_alloc(&AlgebraicNumberType, 0);
  if (self != NULL) {
    if (a) {
      algebraic_number_ops.construct_copy(&self->a, a);
    } else {
      dyadic_rational_t zero;
      dyadic_rational_ops.construct(&zero);
      algebraic_number_ops.construct_from_dyadic_rational(&self->a, &zero);
      dyadic_rational_ops.destruct(&zero);
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
    if (PyUPolynomial_CHECK(f_obj) && PyInt_Check(root_index_obj)) {
      upolynomial_t* f = ((UPolynomialObject*) f_obj)->p;
      long root_index = PyInt_AsLong(root_index_obj);
      size_t roots_count = upolynomial_ops.roots_count(f, 0);
      if (root_index < 0 || root_index >= roots_count) {
        // Not enough roots
        return -1;
      }
      algebraic_number_t* roots = malloc(roots_count * sizeof(algebraic_number_t));
      upolynomial_ops.roots_isolate(f, roots, &roots_count);
      algebraic_number_ops.destruct(&self->a);
      algebraic_number_ops.construct_copy(&self->a, roots + root_index);
      int i;
      for (i = 0; i < roots_count; ++ i) {
        algebraic_number_ops.destruct(roots + i);
      }
      free(roots);
    } else {
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
  algebraic_number_ops.refine(&a->a);
  Py_RETURN_NONE;
}

static int
AlgebraicNumber_cmp(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    // should return -1 and set an exception condition when an error occurred
    return -1;
  }
  // Get arguments
  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;
  // Compare
  return algebraic_number_ops.cmp(&a1->a, &a2->a);
}

static PyObject*
AlgebraicNumber_richcompare(PyObject* self, PyObject* other, int op) {
  PyObject *result = 0;

  if (!PyAlgebraicNumber_CHECK(other)) {
    result = Py_NotImplemented;
  } else {
    algebraic_number_t* self_a = &((AlgebraicNumber*) self)->a;
    algebraic_number_t* other_a = &((AlgebraicNumber*) other)->a;
    int cmp = algebraic_number_ops.cmp(self_a, other_a);

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
  char* cstr = algebraic_number_ops.to_string(&a->a);
  PyObject* pystr = PyString_FromString(cstr);
  free(cstr);
  return pystr;
}

static PyObject*
AlgebraicNumber_add(PyObject* self, PyObject* other) {

  if (PyInt_Check(self)) {
    // To enable summation
    long x = PyInt_AsLong(self);
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

  algebraic_number_t sum;
  algebraic_number_ops.construct_zero(&sum);
  algebraic_number_ops.add(&sum, &a1->a, &a2->a);
  PyObject* result = PyAlgebraicNumber_create(&sum);
  algebraic_number_ops.destruct(&sum);

  return result;
}

static PyObject*
AlgebraicNumber_neg(PyObject* self) {
  return 0;
}

static PyObject*
AlgebraicNumber_sub(PyObject* self, PyObject* other) {
  return 0;
}

static PyObject*
AlgebraicNumber_mul(PyObject* self, PyObject* other) {
  if (!PyAlgebraicNumber_CHECK(self) || !PyAlgebraicNumber_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  AlgebraicNumber* a1 = (AlgebraicNumber*) self;
  AlgebraicNumber* a2 = (AlgebraicNumber*) other;

  algebraic_number_t mul;
  algebraic_number_ops.construct_zero(&mul);
  algebraic_number_ops.mul(&mul, &a1->a, &a2->a);
  PyObject* result = PyAlgebraicNumber_create(&mul);
  algebraic_number_ops.destruct(&mul);

  return result;
}

static PyObject*
AlgebraicNumber_pow(PyObject* self, PyObject* args) {
  return 0;
}
