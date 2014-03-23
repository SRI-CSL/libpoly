/*
 * Assignment.h
 *
 *  Created on: Mar 7, 2014
 *      Author: dejan
 */

#pragma once


#include <Python.h>

#include "assignment.h"

/** Definition of the python Variable object */
typedef struct {
  PyObject_HEAD
  assignment_t* assignment;
} Assignment;

/** Methods on assignments */
extern PyMethodDef Assignment_methods[];

/** Definition of the Assignment type */
extern PyTypeObject AssignmentType;

/** Get the default assignment  */
assignment_t* Assignment_get_default_assignment(void);
