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
%module problems

// Header files for module C++ code
%{
#include "pylith/problems/problemsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator

#include "pylith/problems/Formulation.hh"
#include "pylith/problems/Solver.hh"
#include "pylith/problems/SolverLinear.hh"
#include "pylith/problems/SolverNonlinear.hh"

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

// Interfaces
%include "Formulation.i"
%include "Solver.i"
%include "SolverLinear.i"
%include "SolverNonlinear.i"


// End of file

