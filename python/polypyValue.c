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
#include "utils.h"

#include <structmember.h>
#include <math.h>

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

static long
Value_hash(PyObject* self);

static PyObject*
Value_get_value_between(PyObject* self, PyObject* args);

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

static PyObject*
Value_int(PyObject* self);

static PyObject*
Value_long(PyObject* self);

static PyObject*
Value_float(PyObject* self);

PyMethodDef Value_methods[] = {
    {"to_double", (PyCFunction)Value_to_double, METH_NOARGS, "Returns the approximation of the value"},
    {"get_value_between", (PyCFunction)Value_get_value_between, METH_VARARGS, "Returns a value between this and given value"},
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
     Value_int, // unaryfunc nb_int;
     Value_long, // unaryfunc nb_long;
     Value_float, // unaryfunc nb_float;
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
    PyVarObject_HEAD_INIT(NULL, 0)
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
    &Value_hash,                          /*tp_hash */
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
  Py_TYPE(self)->tp_free((PyObject*)self);
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
    value = INFINITY;
    break;
  case LP_VALUE_MINUS_INFINITY:
    value = -INFINITY;
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
    default:
      assert(0);
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

static long
Value_hash(PyObject* self) {
  Value* v = (Value*) self;
  long hash = lp_value_hash(&v->v);
  if (hash == -1) {
    // value -1 should not be returned as a normal return value
    hash = 0;
  }
  return hash;
}

static PyObject*
Value_get_value_between(PyObject* self, PyObject* args) {
  // self is always a polynomial
  Value* v1 = (Value*) self;
  if (!PyTuple_Check(args) || PyTuple_Size(args) != 1) {
    PyErr_SetString(PyExc_RuntimeError, "get_value_between(): need one argument.");
    return NULL;
  }
  PyObject* other = PyTuple_GetItem(args, 0);

  // other should be a value
  if (!PyValue_CHECK(other)) {
    PyErr_SetString(PyExc_RuntimeError, "get_value_between(): argument not a value.");
    return NULL;
  }
  Value* v2 = (Value*) other;

  // compare the values (they should be different)
  if (v1->v.type == LP_VALUE_NONE || v2->v.type == LP_VALUE_NONE) {
    PyErr_SetString(PyExc_RuntimeError, "get_value_between(): values should not be null.");
    return NULL;
  }
  if (lp_value_cmp(&v1->v, &v2->v) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "get_value_between(): values should be different.");
    return NULL;
  }

  lp_value_t m;
  lp_value_construct_none(&m);
  lp_value_get_value_between(&v1->v, 1, &v2->v, 1, &m);
  PyObject* result = PyValue_create(&m);
  lp_value_destruct(&m);

  return result;
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

// Returns the o converted to an integer object on success, or NULL on failure.
// If the argument is outside the integer range a long object will be returned
// instead. This is the equivalent of the Python expression int(o).
// Return value: New reference.
static PyObject*
Value_int(PyObject* self) {
  Value* value_obj = (Value*) self;
  lp_integer_t int_cast;
  lp_integer_construct(&int_cast);
  lp_value_floor(&value_obj->v, &int_cast);
  PyObject* py_int_cast = integer_to_PyInt(&int_cast);
  lp_integer_destruct(&int_cast);
  return py_int_cast;
}

// Returns the o converted to a long integer object on success, or NULL on
// failure. This is the equivalent of the Python expression long(o).
// Return value: New reference.
static PyObject*
Value_long(PyObject* self) {
  Value* value_obj = (Value*) self;
  lp_integer_t int_cast;
  lp_integer_construct(&int_cast);
  lp_value_floor(&value_obj->v, &int_cast);
  PyObject* py_int_cast = integer_to_PyLong(&int_cast);
  lp_integer_destruct(&int_cast);
  return py_int_cast;
}

// Returns the o converted to a float object on success, or NULL on failure.
// This is the equivalent of the Python expression float(o).
// Return value: New reference.
static PyObject*
Value_float(PyObject* self) {
  Value* value_obj = (Value*) self;
  double value = lp_value_to_double(&value_obj->v);
  return PyFloat_FromDouble(value);
}
