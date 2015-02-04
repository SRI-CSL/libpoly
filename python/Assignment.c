/*
 * Assignment.c
 *
 *  Created on: Mar 7, 2014
 *      Author: dejan
 */

#include "Assignment.h"
#include "Variable.h"
#include "utils.h"

#include <structmember.h>

/** Default variable database */
static lp_assignment_t* default_assignment = 0;

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
    PyObject_HEAD_INIT(NULL)
    0,                            /*ob_size*/
    "polypy.Assignment",          /*tp_name*/
    sizeof(Variable),            /*tp_basicsize*/
    0,                            /*tp_itemsize*/
    (destructor)Assignment_dealloc, /*tp_dealloc*/
    0,                            /*tp_print*/
    0,                            /*tp_getattr*/
    0,                            /*tp_setattr*/
    0,                            /*tp_compare*/
    0,                            /*tp_repr*/
    0,                            /*tp_as_number*/
    0,                            /*tp_as_sequence*/
    0,                            /*tp_as_mapping*/
    0,                            /*tp_hash */
    0,                            /*tp_call*/
    Assignment_str,               /*tp_str*/
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
    Assignment_methods,             /* tp_methods */
    0,                            /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)Assignment_init,      /* tp_init */
    0,                            /* tp_alloc */
    Assignment_new,                 /* tp_new */
};

PyObject*
PyAssignment_create(lp_assignment_t* assignment) {
  Assignment *self;
  self = (Assignment*)AssignmentType.tp_alloc(&AssignmentType, 0);
  if (self != NULL) {
    self->assignment = assignment;
  }
  return (PyObject *)self;
}

static PyObject*
Assignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  return PyAssignment_create(0);
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
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* Assignment_str(PyObject* self) {
  Assignment* a = (Assignment*) self;
  char* a_str = lp_assignment_to_string(a->assignment);
  PyObject* str = PyString_FromString(a_str);
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
      case LP_VALUE_ALGEBRAIC:
        return 0;
        break;
      case LP_VALUE_DYADIC_RATIONAL:
        return dyadic_rational_to_PyFloat(&value->value.dy_q);
        break;
      case LP_VALUE_NONE:
        Py_RETURN_NONE;
        break;
      case LP_VALUE_RATIONAL:
        return 0;
        break;
      default:
        return 0;
      }
    } else {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
}

static PyObject*
Assignment_set_value(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 2) {
      PyObject* var_obj = PyTuple_GetItem(args, 0);
      PyObject* value_obj = PyTuple_GetItem(args, 1);
      if (PyVariable_CHECK(var_obj)) {
        Assignment* a = (Assignment*) self;
        Variable* var = (Variable*) var_obj;
        if (PyFloat_Check(value_obj)) {
          lp_dyadic_rational_t value_dyrat;
          PyFloat_to_dyadic_rational(value_obj, &value_dyrat);
          lp_value_t value;
          lp_value_ops.construct(&value, LP_VALUE_DYADIC_RATIONAL, &value_dyrat);
          lp_assignment_set_value(a->assignment, var->x, &value);
          lp_value_ops.destruct(&value);
          lp_dyadic_rational_destruct(&value_dyrat);
          Py_RETURN_NONE;
        } else if (PyLong_or_Int_Check(value_obj)) {
          lp_integer_t value_int;
          PyLong_or_Int_to_integer(value_obj, lp_Z, &value_int);
          lp_dyadic_rational_t value_dyrat;
          lp_dyadic_rational_construct_from_integer(&value_dyrat, &value_int);
          lp_value_t value;
          lp_value_ops.construct(&value, LP_VALUE_DYADIC_RATIONAL, &value_dyrat);
          lp_assignment_set_value(a->assignment, var->x, &value);
          lp_value_ops.destruct(&value);
          lp_dyadic_rational_destruct(&value_dyrat);
          lp_integer_ops.destruct(&value_int);
          Py_RETURN_NONE;
        } else {
          Py_INCREF(Py_NotImplemented);
          return Py_NotImplemented;
        }
      } else {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
      }
    } else {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
}

static PyObject*
Assignment_unset_value(PyObject* self, PyObject* args) {
  if (PyTuple_Check(args) && PyTuple_Size(args) == 1) {
    PyObject* var_obj = PyTuple_GetItem(args, 0);
    if (PyVariable_CHECK(var_obj)) {
      Variable* var = (Variable*) var_obj;
      Assignment* a = (Assignment*) self;
      lp_assignment_set_value(a->assignment, var->x, 0);
      Py_RETURN_NONE;
    } else {
      Py_INCREF(Py_NotImplemented);
      return Py_NotImplemented;
    }
  } else {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
}
