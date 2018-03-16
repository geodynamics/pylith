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
%module bc

// Header files for module C++ code
%{
#include "pylith/bc/BoundaryCondition.hh"
#include "pylith/bc/Dirichlet.hh"
#include "pylith/bc/DirichletTimeDependent.hh"
#include "pylith/bc/Neumann.hh"
#include "pylith/bc/NeumannTimeDependent.hh"
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
%include "../feassemble/ConstraintPointwise.i"
%include "../feassemble/IntegratorPointwise.i"

%include "BoundaryCondition.i"
%include "Dirichlet.i"
%include "DirichletTimeDependent.i"
%include "Neumann.i"
%include "NeumannTimeDependent.i"



// End of file
