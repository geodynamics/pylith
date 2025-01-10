// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// SWIG interface
%module bc

// Header files for module C++ code
%{
#include "pylith/bc/DirichletTimeDependent.hh"
#include "pylith/bc/NeumannTimeDependent.hh"
#include "pylith/bc/AbsorbingDampers.hh"
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

// Interfaces
%include "../utils/PyreComponent.i"
%include "../problems/Physics.i"

%include "BoundaryCondition.i"
%include "DirichletTimeDependent.i"
%include "NeumannTimeDependent.i"
%include "AbsorbingDampers.i"



// End of file
