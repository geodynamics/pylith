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

// SWIG interface
%module utils

// Header files for module C++ code
%{
#include "pylith/utils/EventLogger.hh"
#include "pylith/utils/TestArray.hh"

#include <petsclog.h> // USES PetscLogEventBegin/End() in inline methods
#include "pylith/utils/arrayfwd.hh" // USES double_array
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

%include "typemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "EventLogger.i"
%include "TestArray.i"

// End of file

