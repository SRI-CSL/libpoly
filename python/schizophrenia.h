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

#include <Python.h>


/**
 *  Returns the char* representation of the python object pyobj.
 *  If pyobj is NULL, then NULL is returned.
 *  If pyobj is a "string" (PyString, PyBytes, or PyUnicode), then the C string 
 *  associated with the object is returned (using UTF-8).
 *  If pyobj is not a string then we try and use the str, then the repr to obtain
 *  a C string version (again using UTF-8 if necessary).
 *  Otherwise we give up and return NULL.
 *
 */
const char* pythonObject2CharStar(PyObject *pyobj);
