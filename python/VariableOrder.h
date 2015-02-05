/*
 * Variable.h
 *
 *  Created on: Feb 12, 2014
 *      Author: dejan
 */

#pragma once

#include <Python.h>

#include <variable_db.h>
#include <variable_order.h>

/** Definition of the python Variable object */
typedef struct {
  PyObject_HEAD
  lp_variable_order_t* var_order;
} VariableOrder;

/** Methods on coefficient rings */
extern PyMethodDef VariableOrder_methods[];

/** Definition of the Variable type */
extern PyTypeObject VariableOrderType;

/** Check if the object is a polynomial */
#define PyVariableOrder_CHECK(arg) \
  (Py_TYPE(arg) == &VariableOrderType)

/**
 * Create a new variable order. Takes over the order reference.
 */
PyObject* VariableOrder_create(lp_variable_order_t* var_order);


/** Set the default variable database */
lp_variable_order_t* VariableOrder_get_default_order(void);

/** Get the variable database */
lp_variable_db_t* Variable_get_default_db(void);
