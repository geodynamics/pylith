// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

// Define the spatialdatacontrib SWIG interface module.

// Set module name
%module materialscontrib

// Header files for module C++ code.
%{
#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "PlaneStrainState.hh"

#include "pylith/utils/arrayfwd.hh"
%}

// Convert standard C++ exceptions to Python exceptions.
%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
} // exception

%include "typemaps.i"
%include "include/doublearray.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "include/numpy.i"
%init %{
import_array();
%}


// Interface files.
%include "materials/Material.i"
%include "materials/ElasticMaterial.i"
%include "PlaneStrainState.i"


// End of file
