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

#include "schizophrenia.h"


const char* pythonObject2CharStar(PyObject *pyobj){
  const char* retval = NULL;
  if(!pyobj){
    return retval;
  } else if(PyBytes_Check(pyobj)){
    retval = PyBytes_AsString(pyobj);
    return retval;
  } else if(PyUnicode_Check(pyobj)) {
    PyObject* str = PyUnicode_AsEncodedString(pyobj, "utf-8", "?");
    retval = PyBytes_AS_STRING(str);
    Py_XDECREF(str);
  } else {
   PyObject* pyob_str = PyObject_Str(pyobj);
   PyObject* str = PyUnicode_AsEncodedString(pyob_str, "utf-8", "?");
   retval = PyBytes_AS_STRING(str);
   Py_XDECREF(pyob_str);
   Py_XDECREF(str);
  }
  return retval;
}
