// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// List of mesh integrators.
%typemap(in) (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >** integrators,
	      const int numIntegrators)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >* integrator = 0;
      int err = SWIG_ConvertPtr(s, (void**) &integrator, 
				$descriptor(pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >*),
				0);
      if (SWIG_IsOK(err))
	$1[i] = (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >*) integrator;
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain mesh integrators.");
	delete[] $1;
	return NULL;
      } // if
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of mesh integrators.");
    return NULL;
  } // if/else
} // typemap(in) [List of mesh integrators.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg) (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::Mesh> >** integrators,
	      const int numIntegrators) {
  delete[] $1;
}

// ----------------------------------------------------------------------
// List of submesh integrators.
%typemap(in) (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >** integrators,
	      const int numIntegrators)
{
  // Check to make sure input is a list.
  if (PyList_Check($input)) {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >*[size] : 0;
    for (int i = 0; i < size; i++) {
      PyObject* s = PyList_GetItem($input,i);
      pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >* integrator = 0;
      int err = SWIG_ConvertPtr(s, (void**) &integrator, 
				$descriptor(pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >*),
				0);
      if (SWIG_IsOK(err))
	$1[i] = (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >*) integrator;
      else {
	PyErr_SetString(PyExc_TypeError, "List must contain submesh integrators.");
	delete[] $1;
	return NULL;
      } // if
    } // for
  } else {
    PyErr_SetString(PyExc_TypeError, "Expected list of submesh integrators.");
    return NULL;
  } // if/else
} // typemap(in) [List of submesh integrators.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg) (pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >** integrators,
	      const int numIntegrators) {
  delete[] $1;
}


// End of file
