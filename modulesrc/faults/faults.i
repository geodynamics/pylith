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
%module faults

// Header files for module C++ code
%{
#include "pylith/faults/SlipTimeFn.hh"
#include "pylith/faults/StepSlipFn.hh"
#include "pylith/faults/ConstRateSlipFn.hh"
#include "pylith/faults/BruneSlipFn.hh"
#include "pylith/faults/LiuCosSlipFn.hh"
#include "pylith/faults/TimeHistorySlipFn.hh"
#include "pylith/faults/EqKinSrc.hh"
#include "pylith/faults/TractPerturbation.hh"
#include "pylith/faults/Fault.hh"
#include "pylith/faults/FaultCohesive.hh"
#include "pylith/faults/FaultCohesiveLagrange.hh"
#include "pylith/faults/FaultCohesiveKin.hh"
#include "pylith/faults/FaultCohesiveDyn.hh"
#include "pylith/faults/FaultCohesiveImpulses.hh"
#include "pylith/faults/FaultCohesiveTract.hh"
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
%include "../feassemble/Integrator.i" // ISA Integrator
%include "../bc/TimeDependent.i" // ISA TimeDependent

%include "SlipTimeFn.i"
%include "StepSlipFn.i"
%include "ConstRateSlipFn.i"
%include "BruneSlipFn.i"
%include "LiuCosSlipFn.i"
%include "TimeHistorySlipFn.i"
%include "EqKinSrc.i"
%include "TractPerturbation.i"
%include "Fault.i"
%include "FaultCohesive.i"
%include "FaultCohesiveLagrange.i"
%include "FaultCohesiveKin.i"
%include "FaultCohesiveDyn.i"
%include "FaultCohesiveImpulses.i"
%include "FaultCohesiveTract.i"

// End of file

