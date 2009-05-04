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
%module bc

// Header files for module C++ code
%{
#include "pylith/bc/BoundaryCondition.hh"
#include "pylith/bc/DirichletBC.hh"
#include "pylith/bc/DirichletBoundary.hh"
#include "pylith/bc/AbsorbingDampers.hh"
#include "pylith/bc/Neumann.hh"
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
%include "../include/doublearray.i"

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
%include "../feassemble/Quadrature.i" // ISA Quadrature
%include "../feassemble/Integrator.i" // ISA Integrator

// template instantiation
%template(SubMeshIntegrator) pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >;

%include "BoundaryCondition.i"
%include "DirichletBC.i"
%include "DirichletBoundary.i"
%include "AbsorbingDampers.i"
%include "Neumann.i"

// End of file

