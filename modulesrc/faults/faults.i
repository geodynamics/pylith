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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module faults

// Header files for module C++ code
%{
#include "pylith/faults/FaultCohesive.hh"
#include "pylith/faults/FaultCohesiveKin.hh"
#include "pylith/faults/KinSrc.hh"
#include "pylith/faults/KinSrcStep.hh"
#include "pylith/faults/KinSrcConstRate.hh"
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
%include "../include/kinsrcarray.i"

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
%include "../feassemble/IntegratorPointwise.i" // ISA IntegratorPointwise

%include "FaultCohesive.i"
%include "FaultCohesiveKin.i"
%include "KinSrc.i"
%include "KinSrcStep.i"
%include "KinSrcConstRate.i"

// End of file

