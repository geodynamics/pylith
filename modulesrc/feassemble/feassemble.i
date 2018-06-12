// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module feassemble

// Header files for module C++ code
%{
#include "pylith/feassemble/IntegratorPointwise.hh"
#include "pylith/feassemble/ConstraintPointwise.hh"
#include "pylith/feassemble/Observer.hh"
#include "pylith/feassemble/ObservedSubject.hh"    
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
}  // exception

%include "typemaps.i"
%include "../include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

%include "../utils/PyreComponent.i"
%include "../topology/FieldBase.i"
%include "Observer.i"
%include "ObservedSubject.i"
%include "IntegratorPointwise.i"
%include "ConstraintPointwise.i"

// End of file
