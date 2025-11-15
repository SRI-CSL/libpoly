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

#pragma once

#include "python.h"

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
