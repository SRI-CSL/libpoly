/*
 * Variable.c
 *
 *  Created on: Feb 12, 2014
 *      Author: dejan
 */

#include "Variable.h"
#include "Integer.h"
#include "Polynomial.h"

#include <structmember.h>

/** Default variable database */
static variable_db_t* default_var_db = 0;

void Variable_init_default_db(void) {
  if (default_var_db) {
    variable_db_ops.detach(default_var_db);
  }
  default_var_db = variable_db_ops.new();
}

variable_db_t* Variable_get_default_db(void) {
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

PyMethodDef Variable_methods[] = {
    {NULL}  /* Sentinel */
};

static PyObject*
Variable_add(PyObject* self, PyObject* args);

static PyObject*
Variable_neg(PyObject* self);

static PyObject*
Variable_sub(PyObject* self, PyObject* args);

static PyObject*
Variable_mul(PyObject* self, PyObject* args);

static PyObject*
Variable_pow(PyObject* self, PyObject* args);

PyNumberMethods Variable_NumberMethods = {
     Variable_add, // binaryfunc nb_add;
     Variable_sub, // binaryfunc nb_subtract;
     Variable_mul, // binaryfunc nb_multiply;
     0, // binaryfunc nb_divide;
     0, // binaryfunc nb_remainder;
     0, // binaryfunc nb_divmod;
     (ternaryfunc)Variable_pow, // ternaryfunc nb_power;
     Variable_neg,   // unaryfunc nb_negative;
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

PyTypeObject VariableType = {
    PyObject_HEAD_INIT(NULL)
    0,                            /*ob_size*/
    "polypy.Variable",            /*tp_name*/
    sizeof(Variable),             /*tp_basicsize*/
    0,                            /*tp_itemsize*/
    (destructor)Variable_dealloc, /*tp_dealloc*/
    0,                            /*tp_print*/
    0,                            /*tp_getattr*/
    0,                            /*tp_setattr*/
    0,                            /*tp_compare*/
    Variable_repr,                /*tp_repr*/
    &Variable_NumberMethods,   /*tp_as_number*/
    0,                            /*tp_as_sequence*/
    0,                            /*tp_as_mapping*/
    0,                            /*tp_hash */
    0,                            /*tp_call*/
    Variable_str,                 /*tp_str*/
    0,                            /*tp_getattro*/
    0,                            /*tp_setattro*/
    0,                            /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES, /*tp_flags*/
    "Variable objects",   /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    Variable_methods,             /* tp_methods */
    0,                            /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)Variable_init,      /* tp_init */
    0,                            /* tp_alloc */
    Variable_new,                 /* tp_new */
};

PyObject*
PyVariable_create(variable_t x) {
  Variable *self;
  self = (Variable*)VariableType.tp_alloc(&VariableType, 0);
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
    if (PyString_Check(obj)) {
      const char* c_str = PyString_AsString(obj);
      variable_t x = variable_db_ops.new_variable(Variable_get_default_db(), c_str);
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
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* Variable_str(PyObject* self) {
  Variable* x = (Variable*) self;
  const char* x_str = variable_db_ops.get_name(Variable_get_default_db(), x->x);
  PyObject* str = PyString_FromString(x_str);
  return str;
}

static PyObject* Variable_repr(PyObject* self) {
  Variable* x = (Variable*) self;
  const char* x_str = variable_db_ops.get_name(Variable_get_default_db(), x->x);
  char* x_repr = malloc(strlen(x_str) + strlen(VariableType.tp_name) + 5);
  sprintf(x_repr, "%s('%s')", VariableType.tp_name, x_str);
  PyObject* str = PyString_FromString(x_repr);
  free(x_repr);
  return str;
}

static
polynomial_t* PyLong_Or_Int_to_polynomial(PyObject* number) {
  const polynomial_context_t* ctx = Polynomial_get_default_context();
  integer_t c;
  PyLong_or_Int_to_integer(number, 0, &c);
  polynomial_t* p_c = polynomial_ops.alloc();
  polynomial_ops.construct_simple(p_c, ctx, &c, 0, 0);
  integer_ops.destruct(&c);
  return p_c;
}

static
polynomial_t* Variable_to_polynomial(PyObject* var) {
  const polynomial_context_t* ctx = Polynomial_get_default_context();
  Variable* x = (Variable*) var;
  integer_t one;
  integer_ops.construct_from_int(Z, &one, 1);
  polynomial_t* p_x = polynomial_ops.alloc();
  polynomial_ops.construct_simple(p_x, ctx, &one, x->x, 1);
  integer_ops.destruct(&one);
  return p_x;
}

static PyObject*
Variable_add_number(PyObject* self, PyObject* other) {

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  polynomial_t* p_sum = polynomial_ops.new(ctx);
  polynomial_ops.add(p_sum, p_x, p_c);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_sum);
}

static PyObject*
Variable_add_Variable(PyObject* self, PyObject* other) {

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_y = Variable_to_polynomial(other);

  // x + c polynomial
  polynomial_t* p_sum = polynomial_ops.new(ctx);
  polynomial_ops.add(p_sum, p_x, p_y);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_y);
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

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);

  // -x polynomial
  polynomial_t* p_neg = polynomial_ops.new(ctx);
  polynomial_ops.neg(p_neg, p_x);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  free(p_x);

  return Polynomial_create(p_neg);
}

static PyObject*
Variable_sub_number(PyObject* self, PyObject* other) {

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  polynomial_t* p_sub = polynomial_ops.new(ctx);
  polynomial_ops.sub(p_sub, p_x, p_c);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_sub);
}

static PyObject*
Variable_sub_Variable(PyObject* self, PyObject* other) {

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_y = Variable_to_polynomial(other);

  // x - y polynomial
  polynomial_t* p_sub = polynomial_ops.new(ctx);
  polynomial_ops.sub(p_sub, p_x, p_y);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_y);
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
    polynomial_ops.neg(result->p, result->p);
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

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_c = PyLong_Or_Int_to_polynomial(other);

  // x + c polynomial
  polynomial_t* p_mul = polynomial_ops.new(ctx);
  polynomial_ops.mul(p_mul, p_x, p_c);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_c);
  free(p_x);
  free(p_c);

  return Polynomial_create(p_mul);
}

static PyObject*
Variable_mul_Variable(PyObject* self, PyObject* other) {

  const polynomial_context_t* ctx = Polynomial_get_default_context();

  // The x polynomial
  polynomial_t* p_x = Variable_to_polynomial(self);
  // The c polynomial
  polynomial_t* p_y = Variable_to_polynomial(other);

  // x + c polynomial
  polynomial_t* p_mul = polynomial_ops.new(ctx);
  polynomial_ops.mul(p_mul, p_x, p_y);

  // Remove temporaries
  polynomial_ops.destruct(p_x);
  polynomial_ops.destruct(p_y);
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
  if (!PyVariable_CHECK(self) || !PyInt_Check(other)) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  } else {
    long n = PyInt_AsLong(other);
    if (n < 0) {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    } else {
      const polynomial_context_t* ctx = Polynomial_get_default_context();
      Variable* var = (Variable*) self;
      integer_t one;
      integer_ops.construct_from_int(Z, &one, 1);
      polynomial_t* pow_x = polynomial_ops.alloc();
      polynomial_ops.construct_simple(pow_x, ctx, &one, var->x, n);
      integer_ops.destruct(&one);
      return Polynomial_create(pow_x);
    }
  }
}
