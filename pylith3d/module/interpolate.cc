// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
//
//  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//
//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.
//
//  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
//  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
//  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include "config.h"

#include <Distribution.hh>
#include <petscmesh.h>
#include <src/dm/mesh/meshpylith.h>
#include <petscmat.h>
#include "journal/debug.h"

#include <Python.h>

#include "interpolate.h"
#include <stdio.h>
#include <string.h>

#include "Numeric/arrayobject.h"

char pypylith3d_interpolatePoints__doc__[] = "";
char pypylith3d_interpolatePoints__name__[] = "interpolatePoints";

PyObject * pypylith3d_interpolatePoints(PyObject *, PyObject *args)
{
  using ALE::Obj;
  PyObject      *pyMesh, *pySol;
  PyArrayObject *pyPoints;

  int ok = PyArg_ParseTuple(args, (char *) "OOO!:interpolatePoints", &pyMesh, &pySol, &PyArray_Type, &pyPoints);
  if (!ok) {
    return 0;
  }
  if ((pyPoints->nd != 2) || (pyPoints->descr->type_num != PyArray_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "points must be a 2d array with double values");
    return 0;
  }
  if (pyPoints->dimensions[1] != 3) {
    PyErr_SetString(PyExc_ValueError, "points must be a 3d");
    return 0;
  }
  if ((pyPoints->strides[0] != 3 * sizeof(double)) || (pyPoints->strides[1] != sizeof(double))) {
    PyErr_SetString(PyExc_ValueError, "points must be a contiguous array");
    return 0;
  }

  Mesh           mesh = (Mesh) PyCObject_AsVoidPtr(pyMesh);
  Vec            sol  = (Vec)  PyCObject_AsVoidPtr(pySol);
  SectionReal    displacement, fullDisplacement;
  const int      numPoints = pyPoints->dimensions[0];
  double        *values;
  PetscErrorCode ierr;

  ierr = MeshGetSectionReal(mesh, "displacement", &displacement);
  ierr = updateDisplacement(displacement, sol);
  ierr = createFullDisplacement(mesh, &fullDisplacement);
  ierr = MeshInterpolatePoints(mesh, fullDisplacement, numPoints, (double *) pyPoints->data, &values);
  ierr = SectionRealDestroy(displacement);
  ierr = PetscFree(values);

  int            dims[2]  = {numPoints, 3};
  PyArrayObject *pyValues = (PyArrayObject *) PyArray_FromDims(2, dims, PyArray_DOUBLE);
  double        *data     = (double *) pyValues->data;

  for(int p = 0; p < numPoints; ++p) {
    for(int d = 0; d < 3; d++) {
      data[p*3+d] = values[p*3+d];
    }
  }

  ierr = PetscFree(values);
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Interpolated points"
    << journal::endl;

  return Py_BuildValue((char *) "N", pyValues);
}

// End of file
