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

#include <petsclog.h> // USES PetscLogEventBegin/End() in inline methods
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

// Interfaces
%include "EventLogger.i"

// End of file

