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

#include "polypyValue.h"

#include <structmember.h>

static void
Value_dealloc(Value* self);

static PyObject*
Value_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
Value_init(Value* self, PyObject* args);

static PyObject*
Value_to_double(PyObject* self);

static int
Value_cmp(PyObject* self, PyObject* args);

static PyObject*
Value_richcompare(PyObject* self, PyObject* other, int op);

static PyObject*
Value_str(PyObject* self);

static PyObject*
Value_add(PyObject* self, PyObject* args);

static PyObject*
Value_neg(PyObject* self);

static PyObject*
Value_sub(PyObject* self, PyObject* args);

static PyObject*
Value_mul(PyObject* self, PyObject* args);

static PyObject*
Value_pow(PyObject* self, PyObject* args);

PyMethodDef Value_methods[] = {
    {"to_double", (PyCFunction)Value_to_double, METH_NOARGS, "Returns the approximation of the value"},
    {NULL}  /* Sentinel */
};

PyNumberMethods Value_NumberMethods = {
     Value_add, // binaryfunc nb_add;
     Value_sub, // binaryfunc nb_subtract;
     Value_mul, // binaryfunc nb_multiply;
     0, // binaryfunc nb_divide;
     0, // binaryfunc nb_remainder;
     0, // binaryfunc nb_divmod;
     (ternaryfunc)Value_pow, // ternaryfunc nb_power;
     Value_neg, // unaryfunc nb_negative;
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

PyTypeObject ValueType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "polypy.Value",   /*tp_name*/
    sizeof(Value),    /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)Value_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    Value_cmp,        /*tp_compare*/
    Value_str,        /*tp_repr*/
    &Value_NumberMethods, /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    Value_str,        /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "Values of different kinds", /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    Value_richcompare,/* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    Value_methods,    /* tp_methods */
    0,                          /* tp_members */
    0,                          /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)Value_init,/* tp_init */
    0,                          /* tp_alloc */
    Value_new,        /* tp_new */
};

static void
Value_dealloc(Value* self)
{
  lp_value_destruct(&self->v);
  self->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyValue_create(const lp_value_t* v) {
  Value *self;
  self = (Value*)ValueType.tp_alloc(&ValueType, 0);
  if (self != NULL) {
    if (v) {
      lp_value_construct_copy(&self->v, v);
    } else {
      lp_value_construct_zero(&self->v);
    }
  }
  return (PyObject *)self;
}

static PyObject*
Value_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyValue_create(0);
}

/** Construct a value from given number. */
static int
Value_init(Value* self, PyObject* args)
{
  assert(0);
  // All fine, initialized
  return 0;
}

static PyObject*
Value_to_double(PyObject* self) {

  Value* a = (Value*) self;

  double value = 0;

  switch (a->v.type) {
  case LP_VALUE_NONE:
    assert(0);
    break;
  case LP_VALUE_INTEGER:
    value = lp_integer_to_int(&a->v.value.z);
    break;
  case LP_VALUE_DYADIC_RATIONAL:
    value = lp_dyadic_rational_to_double(&a->v.value.dy_q);
    break;
  case LP_VALUE_RATIONAL:
    value = lp_rational_to_double(&a->v.value.q);
    break;
  case LP_VALUE_ALGEBRAIC:
    value = lp_algebraic_number_to_double(&a->v.value.a);
    break;
  case LP_VALUE_PLUS_INFINITY:
    assert(0);
    break;
  case LP_VALUE_MINUS_INFINITY:
    assert(0);
    break;
  }

  return PyFloat_FromDouble(value);
}

static int
Value_cmp(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyValue_CHECK(self) || !PyValue_CHECK(other)) {
    // should return -1 and set an exception condition when an error occurred
    return -1;
  }
  // Get arguments
  Value* v1 = (Value*) self;
  Value* v2 = (Value*) other;
  // Compare
  return lp_value_cmp(&v1->v, &v2->v);
}

static PyObject*
Value_richcompare(PyObject* self, PyObject* other, int op) {
  PyObject *result = 0;

  if (!PyValue_CHECK(other)) {
    result = Py_NotImplemented;
  } else {
    lp_value_t* self_v = &((Value*) self)->v;
    lp_value_t* other_v = &((Value*) other)->v;
    int cmp = lp_value_cmp(self_v, other_v);

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

static PyObject* Value_str(PyObject* self) {
  Value* v = (Value*) self;
  char* cstr = lp_value_to_string(&v->v);
  PyObject* pystr = PyString_FromString(cstr);
  free(cstr);
  return pystr;
}

static PyObject*
Value_add(PyObject* self, PyObject* other) {
  assert(0);
  return 0;
}

static PyObject*
Value_neg(PyObject* self) {
  assert(0);
  return 0;
}

static PyObject*
Value_sub(PyObject* self, PyObject* other) {
  assert(0);
  return 0;
}

static PyObject*
Value_mul(PyObject* self, PyObject* other) {
  assert(0);
  return 0;
}

static PyObject*
Value_pow(PyObject* self, PyObject* other) {
  assert(0);
  return 0;
}
