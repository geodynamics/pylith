// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// List of materials.
%typemap(in)(pylith::materials::Material *materials[],
              const int numMaterials)
{
  // Check to make sure input is a list.
  if (PyList_Check($input))
  {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::materials::Material *[size] : 0;
    for (int i = 0; i < size; i++)
    {
      PyObject *s = PyList_GetItem($input, i);
      pylith::materials::Material *material = 0;
      int err = SWIG_ConvertPtr(s, (void **)&material,
                                $descriptor(pylith::materials::Material *),
                                0);
      if (SWIG_IsOK(err))
        $1[i] = (pylith::materials::Material *)material;
      else
      {
        PyErr_SetString(PyExc_TypeError, "List must contain materials.");
        delete[] $1;
        return NULL;
      } // if
    }   // for
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "Expected list of materials.");
    return NULL;
  } // if/else
} // typemap(in) [List of materials.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg)(pylith::materials::Material *materials[],
                   const int numMaterials)
{
  delete[] $1;
}

// ----------------------------------------------------------------------
// List of sources.
%typemap(in)(pylith::sources::Source *sources[],
              const int numSources)
{
  // Check to make sure input is a list.
  if (PyList_Check($input))
  {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::sources::Source *[size] : 0;
    for (int i = 0; i < size; i++)
    {
      PyObject *s = PyList_GetItem($input, i);
      pylith::sources::Source *sources = 0;
      int err = SWIG_ConvertPtr(s, (void **)&sources,
                                $descriptor(pylith::sources::Source *),
                                0);
      if (SWIG_IsOK(err))
        $1[i] = (pylith::sources::Source *)sources;
      else
      {
        PyErr_SetString(PyExc_TypeError, "List must contain sources.");
        delete[] $1;
        return NULL;
      } // if
    }   // for
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "Expected list of sources.");
    return NULL;
  } // if/else
} // typemap(in) [List of sources.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg)(pylith::sources::Source *sources[],
                   const int numSources)
{
  delete[] $1;
}

// ----------------------------------------------------------------------
// List of boundary conditions.
%typemap(in)(pylith::bc::BoundaryCondition *bc[],
              const int numBC)
{
  // Check to make sure input is a list.
  if (PyList_Check($input))
  {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::bc::BoundaryCondition *[size] : 0;
    for (int i = 0; i < size; i++)
    {
      PyObject *s = PyList_GetItem($input, i);
      pylith::bc::BoundaryCondition *bc = 0;
      int err = SWIG_ConvertPtr(s, (void **)&bc,
                                $descriptor(pylith::bc::BoundaryCondition *),
                                0);
      if (SWIG_IsOK(err))
        $1[i] = (pylith::bc::BoundaryCondition *)bc;
      else
      {
        PyErr_SetString(PyExc_TypeError, "List must contain boundary conditions.");
        delete[] $1;
        return NULL;
      } // if
    }   // for
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "Expected list of boundary conditions.");
    return NULL;
  } // if/else
} // typemap(in) [List of bcs.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg)(pylith::bc::BoundaryCondition *bc[],
                   const int numBC)
{
  delete[] $1;
}

// ----------------------------------------------------------------------
// List of initial conditions.
%typemap(in)(pylith::problems::InitialCondition *ic[],
              const int numIC)
{
  // Check to make sure input is a list.
  if (PyList_Check($input))
  {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::problems::InitialCondition *[size] : 0;
    for (int i = 0; i < size; i++)
    {
      PyObject *s = PyList_GetItem($input, i);
      pylith::problems::InitialCondition *ic = 0;
      int err = SWIG_ConvertPtr(s, (void **)&ic,
                                $descriptor(pylith::problems::InitialCondition *),
                                0);
      if (SWIG_IsOK(err))
        $1[i] = (pylith::problems::InitialCondition *)ic;
      else
      {
        PyErr_SetString(PyExc_TypeError, "List must contain boundary conditions.");
        delete[] $1;
        return NULL;
      } // if
    }   // for
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "Expected list of boundary conditions.");
    return NULL;
  } // if/else
} // typemap(in) [List of ics.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg)(pylith::problems::InitialCondition *ic[],
                   const int numIC)
{
  delete[] $1;
}

// ----------------------------------------------------------------------
// List of interfaces.
%typemap(in)(pylith::faults::FaultCohesive *faults[],
              const int numFaults)
{
  // Check to make sure input is a list.
  if (PyList_Check($input))
  {
    const int size = PyList_Size($input);
    $2 = size;
    $1 = (size > 0) ? new pylith::faults::FaultCohesive *[size] : 0;
    for (int i = 0; i < size; i++)
    {
      PyObject *s = PyList_GetItem($input, i);
      pylith::faults::FaultCohesive *fault = 0;
      int err = SWIG_ConvertPtr(s, (void **)&fault,
                                $descriptor(pylith::faults::FaultCohesive *),
                                0);
      if (SWIG_IsOK(err))
        $1[i] = (pylith::faults::FaultCohesive *)fault;
      else
      {
        PyErr_SetString(PyExc_TypeError, "List must contain cohesive faults.");
        delete[] $1;
        return NULL;
      } // if
    }   // for
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "Expected list of cohesive faults.");
    return NULL;
  } // if/else
} // typemap(in) [List of faults.]

// This cleans up the array we malloc'd before the function call
%typemap(freearg)(pylith::faults::FaultCohesive *faults[],
                   const int numFaults)
{
  delete[] $1;
}

// End of file
