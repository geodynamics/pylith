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
%module faults

// Header files for module C++ code
%{
#include "pylith/faults/SlipTimeFn.hh"
#include "pylith/faults/StepSlipFn.hh"
#include "pylith/faults/ConstRateSlipFn.hh"
#include "pylith/faults/BruneSlipFn.hh"
#include "pylith/faults/LiuCosSlipFn.hh"
#include "pylith/faults/EqKinSrc.hh"
#include "pylith/faults/Fault.hh"
#include "pylith/faults/FaultCohesive.hh"
#include "pylith/faults/FaultCohesiveDyn.hh"
#include "pylith/faults/FaultCohesiveKin.hh"

#include "pylith/topology/SubMesh.hh"
#include "pylith/feassemble/Quadrature.hh"
#include "pylith/feassemble/Integrator.hh"
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
%include "../include/chararray.i"
%include "../include/eqkinsrcarray.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "../topology/SubMesh.i" // ISA Integrator<Quadrature<SubMesh> >
%include "../feassemble/Quadrature.i" // ISA Integrator<Quadrature<SubMesh> >
%include "../feassemble/Integrator.i" // ISA Integrator<Quadrature<SubMesh> >

// Template instatiation
%template(SubMeshIntegrator) pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh > >;

%include "SlipTimeFn.i"
%include "StepSlipFn.i"
%include "ConstRateSlipFn.i"
%include "BruneSlipFn.i"
%include "LiuCosSlipFn.i"
%include "EqKinSrc.i"
%include "Fault.i"
%include "FaultCohesive.i"
%include "FaultCohesiveDyn.i"
%include "FaultCohesiveKin.i"

// End of file

