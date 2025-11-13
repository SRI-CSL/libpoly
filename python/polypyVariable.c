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

#include "polypyVariable.h"
#include "polypyInteger.h"
#include "polypyPolynomial.h"
#include "utils.h"

#include <structmember.h>

/** Default variable database */
static lp_variable_db_t* default_var_db = NULL;

void Variable_init_default_db(void) {
  if (default_var_db) {
    lp_variable_db_detach(default_var_db);
  }
  default_var_db = lp_variable_db_new();
}

lp_variable_db_t* Variable_get_default_db(void) {
  if (!default_var_db) {
    Variable_init_default_db();
  }
  return default_var_db;
}

static PyObject*
Variable_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
Variable_init(Variable* self, PyObject* args);

static void
Variable_dealloc(Variable* self);

static PyObject*
Variable_str(PyObject* self);

static PyObject*
Variable_repr(PyObject* self);

static Py_hash_t
Variable_hash(PyObject* self);

static PyObject *
Variable_richcompare(PyObject *self, PyObject *other, int op);

PyMethodDef Variable_methods[] = {
    {NULL}  /* Sentinel */
};

static PyObject*
Variable_add(PyObject* self, PyObject* other);

static PyObject*
Variable_neg(PyObject* self);

static PyObject*
Variable_sub(PyObject* self, PyObject* other);

static PyObject*
Variable_mul(PyObject* self, PyObject* other);

static PyObject*
Variable_pow(PyObject* self, PyObject* other);

PyNumberMethods Variable_NumberMethods = {
     Variable_add,              // binaryfunc nb_add;
     Variable_sub,              // binaryfunc nb_subtract;
     Variable_mul,              // binaryfunc nb_multiply;
     NULL,                      // binaryfunc nb_remainder;
     NULL,                      // binaryfunc nb_divmod;
     (ternaryfunc)Variable_pow, // ternaryfunc nb_power;
     Variable_neg,              // unaryfunc nb_negative;
     NULL,                      // unaryfunc nb_positive;
     NULL,                      // unaryfunc nb_absolute;
     NULL,                      // inquiry nb_bool;
     NULL,                      // unaryfunc nb_invert;
     NULL,                      // binaryfunc nb_lshift;
     NULL,                      // binaryfunc nb_rshift;
     NULL,                      // binaryfunc nb_and;
     NULL,                      // binaryfunc nb_xor;
     NULL,                      // binaryfunc nb_or;
     NULL,                      // unaryfunc nb_int;
     NULL,                      // void *nb_reserved;
     NULL,                      // unaryfunc nb_float;
     NULL,                      // binaryfunc nb_inplace_add;
     NULL,                      // binaryfunc nb_inplace_subtract;
     NULL,                      // binaryfunc nb_inplace_multiply;
     NULL,                      // binaryfunc nb_inplace_remainder;
     NULL,                      // ternaryfunc nb_inplace_power;
     NULL,                      // binaryfunc nb_inplace_lshift;
     NULL,                      // binaryfunc nb_inplace_rshift;
     NULL,                      // binaryfunc nb_inplace_and;
     NULL,                      // binaryfunc nb_inplace_xor;
     NULL,                      // binaryfunc nb_inplace_or;
     NULL,                      // binaryfunc nb_floor_divide;
     NULL,                      // binaryfunc nb_true_divide;
     NULL,                      // binaryfunc nb_inplace_floor_divide;
     NULL,                      // binaryfunc nb_inplace_true_divide;
     NULL,                      // unaryfunc nb_index;
     NULL,                      // binaryfunc nb_matrix_multiply;
     NULL,                      // binaryfunc nb_inplace_matrix_multiply;
};

PyTypeObject VariableType = {
    {PyObject_HEAD_INIT(NULL)},   // PyObject_VAR_HEAD
    "polypy.Variable",            // const char *tp_name;
    sizeof(Variable),             // Py_ssize_t tp_basicsize;
    0,                            // Py_ssize_t tp_itemsize;
    (destructor)Variable_dealloc, // destructor tp_dealloc;
    0,                            // printfunc tp_print;
    NULL,                         // getattrfunc tp_getattr;
    NULL,                         // setattrfunc tp_setattr;
    NULL,                         // PyAsyncMethods *tp_as_async;
    Variable_repr,                // reprfunc tp_repr;
    &Variable_NumberMethods,      // PyNumberMethods *tp_as_number;
    NULL,                         // PySequenceMethods *tp_as_sequence;
    NULL,                         // PyMappingMethods *tp_as_mapping;
    Variable_hash,                // hashfunc tp_hash;
    NULL,                         // ternaryfunc tp_call;
    Variable_str,                 // reprfunc tp_str;
    NULL,                         // getattrofunc tp_getattro;
    NULL,                         // setattrofunc tp_setattro;
    NULL,                         // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,           // unsigned long tp_flags;
    "Variable objects",           // const char *tp_doc;
    NULL,                         // traverseproc tp_traverse;
    NULL,                         // inquiry tp_clear;
    Variable_richcompare,         // richcmpfunc tp_richcompare;
    0,                            // Py_ssize_t tp_weaklistoffset;
    NULL,                         // getiterfunc tp_iter;
    NULL,                         // iternextfunc tp_iternext;
    Variable_methods,             // struct PyMethodDef *tp_methods;
    NULL,                         // istruct PyMemberDef *tp_members;
    NULL,                         // istruct PyGetSetDef *tp_getset;
    NULL,                         // istruct _typeobject *tp_base;
    NULL,                         // iPyObject *tp_dict;
    NULL,                         // idescrgetfunc tp_descr_get;
    NULL,                         // idescrsetfunc tp_descr_set;
    0,                            // iPy_ssize_t tp_dictoffset;
    (initproc)Variable_init,      // initproc tp_init;
    NULL,                         // allocfunc tp_alloc;
    Variable_new,                 // newfunc tp_new;
    NULL,                         // freefunc tp_free;
    NULL,                         // inquiry tp_is_gc;
    NULL,                         // PyObject *tp_bases;
    NULL,                         // PyObject *tp_mro;
    NULL,                         // PyObject *tp_cache;
    NULL,                         // PyObject *tp_subclasses;
    NULL,                         // PyObject *tp_weaklist;
    NULL,                         // destructor tp_del;
    0,                            // unsigned int tp_version_tag;
    NULL,                         // destructor tp_finalize;
};

PyObject*
PyVariable_create(lp_variable_t x) {
  Variable *self = (Variable*)VariableType.tp_alloc(&VariableType, 0);
  if (self != NULL) {
    self->x = x;
  }
  return (PyObject *)self;
}

static PyObject*
Variable_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyVariable_create(0);
}

static int
Variable_init(Variable* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* obj = PyTuple_GetItem(args, 0);
    if (PyUnicode_Check(obj) || PyBytes_Check(obj)) {
      const char* c_str = pythonObject2CharStar(obj);
      lp_variable_t x = lp_variable_db_new_variable(Variable_get_default_db(), c_str);
      self->x = x;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  return 0;
}

static void
Variable_dealloc(Variable* self)
{
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

static PyObject* Variable_str(PyObject* self) {
  Variable* x = (Variable*) self;
  const char* x_str = lp_variable_db_get_name(Variable_get_default_db(), x->x);
  PyObject* str = PyUnicode_FromString(x_str);
  return str;
}

static PyObject* Variable_repr(PyObject* self) {
  Variable* x = (Variable*) self;
  const char* x_str = lp_variable_db_get_name(Variable_get_default_db(), x->x);
  char* x_repr = malloc(strlen(x_str) + strlen(VariableType.tp_name) + 5);
  sprintf(x_repr, "%s('%s')", VariableType.tp_name, x_str);
  PyObject* str = PyUnicode_FromString(x_repr);
  free(x_repr);
  return str;
}

static Py_hash_t Variable_hash(PyObject* self) {
  Variable* x = (Variable*) self;
  return x->x;
}

static
lp_polynomial_t* PyLong_Or_Int_to_polynomial(PyObject* number) {
  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();
  lp_integer_t c;
  PyLong_or_Int_to_integer(number, NULL, &c);
  lp_polynomial_t* p_c = lp_polynomial_alloc();
  lp_polynomial_construct_simple(p_c, ctx, &c, 0, 0);
  lp_integer_destruct(&c);
  return p_c;
}

static
lp_polynomial_t* Variable_to_polynomial(PyObject* var) {
  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();
  Variable* x = (Variable*) var;
  lp_integer_t one;
  lp_integer_construct_from_int(lp_Z, &one, 1);
  lp_polynomial_t* p_x = lp_polynomial_alloc();
  lp_polynomial_construct_simple(p_x, ctx, &one, x->x, 1);
  lp_integer_destruct(&one);
  return p_x;
}

static PyObject*
Variable_add_number(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  lp_polynomial_t* p_sum = lp_polynomial_new(ctx);
  lp_polynomial_add(p_sum, p_x, p_c);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_sum);
}

static PyObject*
Variable_add_Variable(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_y = Variable_to_polynomial(other);

  // x + c polynomial
  lp_polynomial_t* p_sum = lp_polynomial_new(ctx);
  lp_polynomial_add(p_sum, p_x, p_y);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_y);
  free(p_x);
  free(p_y);

  return Polynomial_create(p_sum);
}

static PyObject*
Variable_add(PyObject* self, PyObject* other) {
  // Integer addition
  if (PyLong_or_Int_Check(other)) {
    return Variable_add_number(self, other);
  } else if (PyLong_or_Int_Check(self)) {
    return Variable_add_number(other, self);
  } else if (PyVariable_CHECK(other)) {
    return Variable_add_Variable(self, other);
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
}

static PyObject*
Variable_neg(PyObject* self) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);

  // -x polynomial
  lp_polynomial_t* p_neg = lp_polynomial_new(ctx);
  lp_polynomial_neg(p_neg, p_x);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  free(p_x);

  return Polynomial_create(p_neg);
}

static PyObject*
Variable_sub_number(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  lp_polynomial_t* p_sub = lp_polynomial_new(ctx);
  lp_polynomial_sub(p_sub, p_x, p_c);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_sub);
}

static PyObject*
Variable_sub_Variable(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_y = Variable_to_polynomial(other);

  // x - y polynomial
  lp_polynomial_t* p_sub = lp_polynomial_new(ctx);
  lp_polynomial_sub(p_sub, p_x, p_y);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_y);
  free(p_x);
  free(p_y);

  return Polynomial_create(p_sub);
}

static PyObject*
Variable_sub(PyObject* self, PyObject* other) {
  // Integer addition
  if (PyLong_or_Int_Check(other)) {
    return Variable_sub_number(self, other);
  } else if (PyLong_or_Int_Check(self)) {
    Polynomial* result = (Polynomial*) Variable_sub_number(other, self);
    lp_polynomial_neg(result->p, result->p);
    return (PyObject*) result;
  } else if (PyVariable_CHECK(other)) {
    return Variable_sub_Variable(self, other);
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
}

static PyObject*
Variable_mul_number(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  lp_polynomial_t* p_mul = lp_polynomial_new(ctx);
  lp_polynomial_mul(p_mul, p_x, p_c);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_mul);
}

static PyObject*
Variable_mul_Variable(PyObject* self, PyObject* other) {

  const lp_polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  lp_polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  lp_polynomial_t* p_y = Variable_to_polynomial(other);

  // x + c polynomial
  lp_polynomial_t* p_mul = lp_polynomial_new(ctx);
  lp_polynomial_mul(p_mul, p_x, p_y);

  // Remove temporaries
  lp_polynomial_destruct(p_x);
  lp_polynomial_destruct(p_y);
  free(p_x);
  free(p_y);

  return Polynomial_create(p_mul);
}

static PyObject*
Variable_mul(PyObject* self, PyObject* other) {
  // Integer addition
  if (PyLong_or_Int_Check(other)) {
    return Variable_mul_number(self, other);
  } else if (PyLong_or_Int_Check(self)) {
    return Variable_mul_number(other, self);
  } else if (PyVariable_CHECK(other)) {
    return Variable_mul_Variable(other, self);
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
}

static PyObject*
Variable_pow(PyObject* self, PyObject* other) {
  // Check arguments
  if (!PyVariable_CHECK(self) || !PyLong_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  } else {
    long n = PyLong_AsLong(other);
    if (n < 0) {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    } else {
      const lp_polynomial_context_t* ctx = Polynomial_get_default_context();
      Variable* var = (Variable*) self;
      lp_integer_t one;
      lp_integer_construct_from_int(lp_Z, &one, 1);
      lp_polynomial_t* pow_x = lp_polynomial_alloc();
      lp_polynomial_construct_simple(pow_x, ctx, &one, var->x, n);
      lp_integer_destruct(&one);
      return Polynomial_create(pow_x);
    }
  }
}

static PyObject *
Variable_richcompare(PyObject *self, PyObject *other, int op){
  PyObject *result = Py_NotImplemented;

  if (!PyVariable_CHECK(other)) {
    if(op == Py_EQ){ return Py_False; }
    if(op == Py_NE){ return Py_True; }
  } else {
    Variable* x = (Variable*) self;
    Variable* y = (Variable*) other;

    switch (op) {
    case Py_LT:
      result = (x->x < y->x ? Py_True : Py_False);
      break;
    case Py_LE:
      result = (x->x <= y->x ? Py_True : Py_False);
      break;
    case Py_EQ:
      result = (x->x == y->x ? Py_True : Py_False);
      break;
    case Py_NE:
      result = (x->x != y->x ? Py_True : Py_False);
      break;
    case Py_GT:
      result = (x->x > y->x ? Py_True : Py_False);
      break;
    case Py_GE:
      result = (x->x >= y->x ? Py_True : Py_False);
      break;
    default:
      assert(0);
    }
  }
  return result;
}

