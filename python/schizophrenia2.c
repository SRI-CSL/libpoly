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
  } else if(PyString_Check(pyobj)){
    retval = PyString_AsString(pyobj);
    return retval;
  } else {
    PyObject* pyobj_str = PyObject_Str(pyobj);
    retval = PyString_AsString(pyobj_str);
  }
  return retval;
}

