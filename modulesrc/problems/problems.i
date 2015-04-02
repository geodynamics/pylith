// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module problems

// Header files for module C++ code
%{
#include "pylith/problems/problemsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator

#include "pylith/problems/Formulation.hh"
#include "pylith/problems/Explicit.hh"
#include "pylith/problems/Implicit.hh"
#include "pylith/problems/Solver.hh"
#include "pylith/problems/SolverLinear.hh"
#include "pylith/problems/SolverNonlinear.hh"
#include "pylith/problems/SolverLumped.hh"
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
%include "../include/integratorarray.i"
%include "../include/scalartypemaps.i"

// Interfaces
%include "Formulation.i"
%include "Explicit.i"
%include "Implicit.i"
%include "Solver.i"
%include "SolverLinear.i"
%include "SolverNonlinear.i"
%include "SolverLumped.i"


// End of file

