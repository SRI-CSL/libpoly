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

#include "polypyAssignment.h"
#include "polypyVariable.h"
#include "utils.h"
#include "polypyAlgebraicNumber.h"
#include "polypyValue.h"

/** Default variable database */
static lp_assignment_t* default_assignment = NULL;

lp_assignment_t* Assignment_get_default_assignment(void) {
  if (!default_assignment) {
    default_assignment = lp_assignment_new(Variable_get_default_db());
  }
  return default_assignment;
}

static PyObject*
Assignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
Assignment_init(Assignment* self, PyObject* args);

static void
Assignment_dealloc(Assignment* self);

static PyObject*
Assignment_str(PyObject* self);

static PyObject*
Assignment_get_value(PyObject* self, PyObject* args);

static PyObject*
Assignment_set_value(PyObject* self, PyObject* args);

static PyObject*
Assignment_unset_value(PyObject* self, PyObject* args);

PyMethodDef Assignment_methods[] = {
    {"get_value", (PyCFunction)Assignment_get_value, METH_VARARGS, "Returns the value of the variable"},
    {"set_value", (PyCFunction)Assignment_set_value, METH_VARARGS, "Sets the value of the variable"},
    {"unset_value", (PyCFunction)Assignment_unset_value, METH_VARARGS, "Marks the variable as unassigned"},
    {NULL}  /* Sentinel */
};


PyTypeObject AssignmentType = {
    {PyObject_HEAD_INIT(NULL)},    // PyObject_VAR_HEAD
    "polypy.Assignment",            // const char *tp_name;
    sizeof(Variable),               // Py_ssize_t tp_basicsize; 
    0,                              // Py_ssize_t tp_itemsize;
    (destructor)Assignment_dealloc, // destructor tp_dealloc;
    0,                              // printfunc tp_print;
    NULL,                           // getattrfunc tp_getattr;
    NULL,                           // setattrfunc tp_setattr;
    NULL,                           // PyAsyncMethods *tp_as_async;
    NULL,                           // reprfunc tp_repr;
    NULL,                           // PyNumberMethods *tp_as_number;
    NULL,                           // PySequenceMethods *tp_as_sequence;
    NULL,                           // PyMappingMethods *tp_as_mapping;
    NULL,                           // hashfunc tp_hash;
    NULL,                           // ternaryfunc tp_call;
    Assignment_str,                 // reprfunc tp_str;
    NULL,                           // getattrofunc tp_getattro;
    NULL,                           // setattrofunc tp_setattro;
    NULL,                           // PyBufferProcs *tp_as_buffer;
    Py_TPFLAGS_DEFAULT,             // unsigned long tp_flags;
    "Variable objects",             // const char *tp_doc;
    NULL,                           // traverseproc tp_traverse;
    NULL,                           // inquiry tp_clear;
    NULL,                           // richcmpfunc tp_richcompare;
    0,                              // Py_ssize_t tp_weaklistoffset;
    NULL,                           // getiterfunc tp_iter;
    NULL,                           // iternextfunc tp_iternext;
    Assignment_methods,             // struct PyMethodDef *tp_methods;
    NULL,                           // struct PyMemberDef *tp_members;
    NULL,                           // struct PyGetSetDef *tp_getset;
    NULL,                           // struct _typeobject *tp_base;
    NULL,                           // PyObject *tp_dict;
    NULL,                           // descrgetfunc tp_descr_get;
    NULL,                           // descrsetfunc tp_descr_set;
    0,                              // Py_ssize_t tp_dictoffset;
    (initproc)Assignment_init,      // initproc tp_init;
    NULL,                           // allocfunc tp_alloc;
    Assignment_new,                 // newfunc tp_new;
    NULL,                           // freefunc tp_free;
    NULL,                           // inquiry tp_is_gc;
    NULL,                           // PyObject *tp_bases;
    NULL,                           // PyObject *tp_mro;
    NULL,                           // PyObject *tp_cache;
    NULL,                           // PyObject *tp_subclasses;
    NULL,                           // PyObject *tp_weaklist;
    NULL,                           // destructor tp_del;
    0,                              // unsigned int tp_version_tag;
    NULL,                           // destructor tp_finalize;
};

PyObject*
PyAssignment_create(lp_assignment_t* assignment) {
  Assignment *self = (Assignment*)AssignmentType.tp_alloc(&AssignmentType, 0);
  if (self != NULL) {
    self->assignment = assignment;
  }
  return (PyObject *)self;
}

static PyObject*
Assignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyAssignment_create(NULL);
}

static int
Assignment_init(Assignment* self, PyObject* args)
{
  if (PyTuple_Check(args) && PyTuple_Size(args) == 0) {
    self->assignment = lp_assignment_new(Variable_get_default_db());
  } else {
    return -1;
  }
  return 0;
}

static void
Assignment_dealloc(Assignment* self)
{
  if (self->assignment) {
    lp_assignment_delete(self->assignment);
  }
  ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

static PyObject* Assignment_str(PyObject* self) {
  Assignment* a = (Assignment*) self;
  char* a_str = lp_assignment_to_string(a->assignment);
  PyObject* str = PyUnicode_FromString(a_str);
  free(a_str);
  return str;
}

static PyObject*
Assignment_get_value(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* var_obj = PyTuple_GetItem(args, 0);
    if (PyVariable_CHECK(var_obj)) {
      Variable* var = (Variable*) var_obj;
      Assignment* a = (Assignment*) self;
      const lp_value_t* value = lp_assignment_get_value(a->assignment, var->x);
      switch (value->type) {
      case LP_VALUE_DYADIC_RATIONAL:
        return dyadic_rational_to_PyFloat(&value->value.dy_q);
      case LP_VALUE_NONE:
        Py_RETURN_NONE;
      case LP_VALUE_RATIONAL:
      case LP_VALUE_ALGEBRAIC:
      default:
        return NULL;
      }
    }
  }
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}

static PyObject*
Assignment_set_value(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args)) {
    if (PyTuple_Size(args) == 2) {
      PyObject* var_obj = PyTuple_GetItem(args, 0);
      PyObject* value_obj = PyTuple_GetItem(args, 1);
      if (PyVariable_CHECK(var_obj)) {
        Assignment* a = (Assignment*) self;
        Variable* var = (Variable*) var_obj;
        if (PyFloat_Check(value_obj)) {
          lp_dyadic_rational_t value_dyrat;
          PyFloat_to_dyadic_rational(value_obj, &value_dyrat);
          lp_value_t value;
          lp_value_construct(&value, LP_VALUE_DYADIC_RATIONAL, &value_dyrat);
          lp_assignment_set_value(a->assignment, var->x, &value);
          lp_value_destruct(&value);
          lp_dyadic_rational_destruct(&value_dyrat);
          Py_RETURN_NONE;
        } else if (PyLong_or_Int_Check(value_obj)) {
          lp_integer_t value_int;
          PyLong_or_Int_to_integer(value_obj, lp_Z, &value_int);
          lp_dyadic_rational_t value_dyrat;
          lp_dyadic_rational_construct_from_integer(&value_dyrat, &value_int);
          lp_value_t value;
          lp_value_construct(&value, LP_VALUE_DYADIC_RATIONAL, &value_dyrat);
          lp_assignment_set_value(a->assignment, var->x, &value);
          lp_value_destruct(&value);
          lp_dyadic_rational_destruct(&value_dyrat);
          lp_integer_destruct(&value_int);
          Py_RETURN_NONE;
        } else if (PyAlgebraicNumber_CHECK(value_obj)) {
          lp_value_t value;
          lp_value_construct(&value, LP_VALUE_ALGEBRAIC, &((AlgebraicNumber*) value_obj)->a);
          lp_assignment_set_value(a->assignment, var->x, &value);
          lp_value_destruct(&value);
          Py_RETURN_NONE;
        } else if (PyValue_CHECK(value_obj)) {
          Value* value = (Value*) value_obj;
          lp_assignment_set_value(a->assignment, var->x, &value->v);
          Py_RETURN_NONE;
        } else {
          PyErr_SetString(PyExc_RuntimeError, "set_value(): not a value.");
          return NULL;
        }
      } else {
        PyErr_SetString(PyExc_RuntimeError, "set_value(): not a variable.");
        return NULL;
      }
    } else if (PyTuple_Size(args) == 3) {
      // Three arguments, rational
      PyObject* var_obj = PyTuple_GetItem(args, 0);
      PyObject* p_obj = PyTuple_GetItem(args, 1);
      PyObject* q_obj = PyTuple_GetItem(args, 2);
      if (PyLong_or_Int_Check(p_obj) && PyLong_or_Int_Check(q_obj)) {
        Assignment* a = (Assignment*) self;
        Variable* var = (Variable*) var_obj;
        lp_integer_t p_int, q_int;
        PyLong_or_Int_to_integer(p_obj, lp_Z, &p_int);
        PyLong_or_Int_to_integer(q_obj, lp_Z, &q_int);
        lp_rational_t value_rat;
        lp_rational_construct_from_div(&value_rat, &p_int, &q_int);
        lp_value_t value;
        lp_value_construct(&value, LP_VALUE_RATIONAL, &value_rat);
        lp_assignment_set_value(a->assignment, var->x, &value);
        lp_value_destruct(&value);
        lp_rational_destruct(&value_rat);
        lp_integer_destruct(&p_int);
        lp_integer_destruct(&q_int);
        Py_RETURN_NONE;
      } else {
        PyErr_SetString(PyExc_RuntimeError, "set_value(): for rationals, both numerator and denominator must be integer.");
        return NULL;
      }
    } else {
      PyErr_SetString(PyExc_RuntimeError, "set_value(): need two or three arguments.");
      return NULL;
    }
  } else {
    PyErr_SetString(PyExc_RuntimeError, "set_value(): need two or three arguments.");
    return NULL;
  }
}

static PyObject*
Assignment_unset_value(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* var_obj = PyTuple_GetItem(args, 0);
    if (PyVariable_CHECK(var_obj)) {
      Variable* var = (Variable*) var_obj;
      Assignment* a = (Assignment*) self;
      lp_assignment_set_value(a->assignment, var->x, NULL);
      Py_RETURN_NONE;
    } else {
      PyErr_SetString(PyExc_RuntimeError, "set_value(): not a variable.");
      return NULL;
    }
  } else {
    PyErr_SetString(PyExc_RuntimeError, "set_value(): need one argument.");
    return NULL;
  }
}
