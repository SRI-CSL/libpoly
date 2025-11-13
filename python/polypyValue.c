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

#include "polypyValue.h"
#include "utils.h"
#include "polypyAlgebraicNumber.h"

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

static PyObject*
Value_richcompare(PyObject* self, PyObject* other, int op);

static PyObject*
Value_str(PyObject* self);

static Py_hash_t
Value_hash(PyObject* self);

static PyObject*
Value_get_value_between(PyObject* self, PyObject* args);

static PyObject*
Value_add(PyObject* self, PyObject* other);

static PyObject*
Value_neg(PyObject* self);

static PyObject*
Value_sub(PyObject* self, PyObject* other);

static PyObject*
Value_mul(PyObject* self, PyObject* other);

static PyObject*
Value_div(PyObject* self, PyObject* other);

static PyObject*
Value_pow(PyObject* self, PyObject* other);

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
     Value_add,              // binaryfunc nb_add;
     Value_sub,              // binaryfunc nb_subtract;
     Value_mul,              // binaryfunc nb_multiply;
     NULL,                   // binaryfunc nb_remainder;
     NULL,                   // binaryfunc nb_divmod;
     (ternaryfunc)Value_pow, // ternaryfunc nb_power;
     Value_neg,              // unaryfunc nb_negative;
     NULL,                   // unaryfunc nb_positive;
     NULL,                   // unaryfunc nb_absolute;
     NULL,                   // inquiry nb_bool;
     NULL,                   // unaryfunc nb_invert;
     NULL,                   // binaryfunc nb_lshift;
     NULL,                   // binaryfunc nb_rshift;
     NULL,                   // binaryfunc nb_and;
     NULL,                   // binaryfunc nb_xor;
     NULL,                   // binaryfunc nb_or;
     Value_long,             // unaryfunc nb_int;
     NULL,                   // void *nb_reserved;
     Value_float,            // unaryfunc nb_float;
     NULL,                   // binaryfunc nb_inplace_add;
     NULL,                   // binaryfunc nb_inplace_subtract;
     NULL,                   // binaryfunc nb_inplace_multiply;
     NULL,                   // binaryfunc nb_inplace_remainder;
     NULL,                   // ternaryfunc nb_inplace_power;
     NULL,                   // binaryfunc nb_inplace_lshift;
     NULL,                   // binaryfunc nb_inplace_rshift;
     NULL,                   // binaryfunc nb_inplace_and;
     NULL,                   // binaryfunc nb_inplace_xor;
     NULL,                   // binaryfunc nb_inplace_or;
     NULL,                   // binaryfunc nb_floor_divide;
     Value_div,              // binaryfunc nb_true_divide;
     NULL,                   // binaryfunc nb_inplace_floor_divide;
     NULL,                   // binaryfunc nb_inplace_true_divide;
     NULL,                   // unaryfunc nb_index;
     NULL,                   // binaryfunc nb_matrix_multiply;
     NULL,                   // binaryfunc nb_inplace_matrix_multiply;
};

PyTypeObject ValueType = {
    {PyObject_HEAD_INIT(NULL)},   // PyObject_VAR_HEAD
    "polypy.Value",              // const char *tp_name;
    sizeof(Value),               // Py_ssize_t tp_basicsize;
    0,                           // Py_ssize_t tp_itemsize;
    (destructor)Value_dealloc,   // destructor tp_dealloc;
    0,                           // printfunc tp_print;
    NULL,                        // getattrfunc tp_getattr;
    NULL,                        // setattrfunc tp_setattr;
    NULL,                        // PyAsyncMethods *tp_as_async;
    Value_str,                   // reprfunc tp_repr;
    &Value_NumberMethods,        // PyNumberMethods *tp_as_number;
    NULL,                        // PySequenceMethods *tp_as_sequence;
    NULL,                        // PyMappingMethods *tp_as_mapping;
    &Value_hash,                 // hashfunc tp_hash;
    NULL,                        // ternaryfunc tp_call;
    Value_str,                   // reprfunc tp_str;
    NULL,                        // getattrofunc tp_getattro;
    NULL,                        // setattrofunc tp_setattro;
    NULL,                        // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,          // unsigned long tp_flags;
    "Values of different kinds", // const char *tp_doc;
    NULL,                        // traverseproc tp_traverse;
    NULL,                        // inquiry tp_clear;
    Value_richcompare,           // richcmpfunc tp_richcompare;
    0,                           // Py_ssize_t tp_weaklistoffset;
    NULL,                        // getiterfunc tp_iter;
    NULL,                        // iternextfunc tp_iternext;
    Value_methods,               // struct PyMethodDef *tp_methods;
    NULL,                        // struct PyMemberDef *tp_members;
    NULL,                        // struct PyGetSetDef *tp_getset;
    NULL,                        // struct _typeobject *tp_base;
    NULL,                        // PyObject *tp_dict;
    NULL,                        // descrgetfunc tp_descr_get;
    NULL,                        // descrsetfunc tp_descr_set;
    0,                           // Py_ssize_t tp_dictoffset;
    (initproc)Value_init,        // initproc tp_init;
    NULL,                        // allocfunc tp_alloc;
    Value_new,                   // newfunc tp_new;
    NULL,                        // freefunc tp_free;
    NULL,                        // inquiry tp_is_gc;
    NULL,                        // PyObject *tp_bases;
    NULL,                        // PyObject *tp_mro;
    NULL,                        // PyObject *tp_cache;
    NULL,                        // PyObject *tp_subclasses;
    NULL,                        // PyObject *tp_weaklist;
    NULL,                        // destructor tp_del;
    0,                           // unsigned int tp_version_tag;
    NULL,                        // destructor tp_finalize;
};

static void
Value_dealloc(Value* self)
{
  lp_value_destruct(&self->v);
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

PyObject*
PyValue_create(const lp_value_t* v) {
  Value *self = (Value*)ValueType.tp_alloc(&ValueType, 0);
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
  return PyValue_create(NULL);
}

/** Construct a value from given number. */
static int
Value_init(Value* self, PyObject* args)
{
  if (PyTuple_Check(args)) {
    if (PyTuple_Size(args) == 0) {
      lp_value_construct_zero(&self->v);
    } else if (PyTuple_Size(args) == 1) {
      PyObject* v = PyTuple_GetItem(args, 0);
      if (PyLong_Check(v)) {
        long v_int = PyLong_AsLong(v);
        lp_value_construct_int(&self->v, v_int);
      } else if (PyAlgebraicNumber_CHECK(v)) {
        AlgebraicNumber* v_alg = (AlgebraicNumber*) v;
        lp_value_construct(&self->v, LP_VALUE_ALGEBRAIC, &v_alg->a);
      } else {
        return -1;
      }
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

static PyObject*
Value_richcompare(PyObject* self, PyObject* other, int op) {
  PyObject *result = NULL;

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
  PyObject* pystr = PyUnicode_FromString(cstr);
  free(cstr);
  return pystr;
}

static Py_hash_t
Value_hash(PyObject* self) {
  Value* v = (Value*) self;
  Py_hash_t hash = lp_value_hash(&v->v);
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
  if (!PyValue_CHECK(self) || !PyValue_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v1 = (Value*) self;
  Value* v2 = (Value*) other;

  lp_value_t add;
  lp_value_construct_none(&add);
  lp_value_add(&add, &v1->v, &v2->v);
  PyObject* result = PyValue_create(&add);
  lp_value_destruct(&add);

  return result;
}

static PyObject*
Value_neg(PyObject* self) {
  if (!PyValue_CHECK(self)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v = (Value*) self;

  lp_value_t neg;
  lp_value_construct_none(&neg);
  lp_value_neg(&neg, &v->v);
  PyObject* result = PyValue_create(&neg);
  lp_value_destruct(&neg);

  return result;
}

static PyObject*
Value_sub(PyObject* self, PyObject* other) {
  if (!PyValue_CHECK(self) || !PyValue_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v1 = (Value*) self;
  Value* v2 = (Value*) other;

  lp_value_t sub;
  lp_value_construct_none(&sub);
  lp_value_sub(&sub, &v1->v, &v2->v);
  PyObject* result = PyValue_create(&sub);
  lp_value_destruct(&sub);

  return result;
}

static PyObject*
Value_mul(PyObject* self, PyObject* other) {
  if (!PyValue_CHECK(self) || !PyValue_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v1 = (Value*) self;
  Value* v2 = (Value*) other;

  lp_value_t mul;
  lp_value_construct_none(&mul);
  lp_value_mul(&mul, &v1->v, &v2->v);
  PyObject* result = PyValue_create(&mul);
  lp_value_destruct(&mul);

  return result;
}

static PyObject*
Value_div(PyObject* self, PyObject* other) {
  if (!PyValue_CHECK(self) || !PyValue_CHECK(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v1 = (Value*) self;
  Value* v2 = (Value*) other;

  lp_value_t div;
  lp_value_construct_none(&div);
  lp_value_div(&div, &v1->v, &v2->v);
  PyObject* result = PyValue_create(&div);
  lp_value_destruct(&div);

  return result;
}


static PyObject*
Value_pow(PyObject* self, PyObject* other) {
  if (!PyValue_CHECK(self) || !PyLong_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }

  Value* v = (Value*) self;
  long n = PyLong_AsLong(other);

  lp_value_t mul;
  lp_value_construct_none(&mul);
  lp_value_pow(&mul, &v->v, n);
  PyObject* result = PyValue_create(&mul);
  lp_value_destruct(&mul);

  return result;
}

// Returns the o converted to a long integer object on success, or NULL on
// failure. This is the equivalent of the Python 3 expression int(o).
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
