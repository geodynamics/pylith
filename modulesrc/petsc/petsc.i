// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
// SWIG interface
%module petsc

// Header files for module C++ code
%{
#include "pylith/petsc/Application.hh"
#include "pylith/petsc/EventLogger.hh"
#include "pylith/petsc/Options.hh"
#include "pylith/petsc/PetscVersion.hh"

#include <petsclog.h> // USES PetscLogEventBegin/End() in inline methods
%}

%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& err) {
        SWIG_exception (SWIG_RuntimeError, err.what ());
    } // try/catch
} // exception

%include "typemaps.i"
%include "../include/scalartypemaps.i"
%include "../include/chararray.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
    import_array();
%}

// Interfaces
%include "PetscVersion.i"
%include "Application.i"
%include "Options.i"
%include "EventLogger.i"

// End of file
