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

// SWIG interface
%module utils

// Header files for module C++ code
%{
#include "pylith/utils/EventLogger.hh"
#include "pylith/utils/PyreComponent.hh"
#include "pylith/utils/PylithVersion.hh"
#include "pylith/utils/PetscVersion.hh"
#include "pylith/utils/DependenciesVersion.hh"
#include "pylith/utils/TestArray.hh"
#include "pylith/utils/constdefs.h"

#include <petsclog.h> // USES PetscLogEventBegin/End() in inline methods
#include "pylith/utils/arrayfwd.hh" // USES scalar_array
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
%include "../include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "pylith_general.i"
%include "EventLogger.i"
%include "PyreComponent.i"
%include "PylithVersion.i"
%include "PetscVersion.i"
%include "DependenciesVersion.i"
%include "TestArray.i"
%include "constdefs.i"

// End of file
