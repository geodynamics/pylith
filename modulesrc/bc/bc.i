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
%module bc

// Header files for module C++ code
%{
#include "pylith/bc/BoundaryCondition.hh"
#include "pylith/bc/BoundaryConditionPoints.hh"
#include "pylith/bc/BCIntegratorSubMesh.hh"
#include "pylith/bc/TimeDependent.hh"
#include "pylith/bc/TimeDependentPoints.hh"
#include "pylith/bc/DirichletBC.hh"
#include "pylith/bc/DirichletBoundary.hh"
#include "pylith/bc/AbsorbingDampers.hh"
#include "pylith/bc/Neumann.hh"
#include "pylith/bc/PointForce.hh"
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
%include "../feassemble/Constraint.i" // ISA Constraint
%include "../feassemble/Integrator.i" // ISA Integrator

%include "BoundaryCondition.i"
%include "BoundaryConditionPoints.i"
%include "BCIntegratorSubMesh.i"
%include "TimeDependent.i"
%include "TimeDependentPoints.i"
%include "DirichletBC.i"
%include "DirichletBoundary.i"
%include "AbsorbingDampers.i"
%include "Neumann.i"
%include "PointForce.i"

// End of file

