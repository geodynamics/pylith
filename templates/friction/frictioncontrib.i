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

// Define the frictioncontrib SWIG interface module.

// Set module name
%module frictioncontrib

// Header files for module C++ code.
%{
#include "pylith/friction/frictionfwd.hh" // forward declarations

#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "ViscousFriction.hh"

#include "pylith/utils/types.hh"
#include "pylith/utils/array.hh"
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
%include "include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "include/numpy.i"
%init %{
import_array();
%}


// Interface files.
%include "friction/FrictionModel.i"
%include "ViscousFriction.i"


// End of file
