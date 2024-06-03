// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// SWIG interface
%module faults

// Header files for module C++ code
%{
#include "pylith/faults/FaultCohesive.hh"
#include "pylith/faults/FaultCohesiveKin.hh"
#include "pylith/faults/FaultCohesiveImpulses.hh"
#include "pylith/faults/KinSrc.hh"
#include "pylith/faults/KinSrcStep.hh"
#include "pylith/faults/KinSrcRamp.hh"
#include "pylith/faults/KinSrcConstRate.hh"
#include "pylith/faults/KinSrcBrune.hh"
#include "pylith/faults/KinSrcLiuCos.hh"
#include "pylith/faults/KinSrcTimeHistory.hh"
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
%include "../problems/Physics.i"

%include "FaultCohesive.i"
%include "FaultCohesiveKin.i"
%include "FaultCohesiveImpulses.i"
%include "KinSrc.i"
%include "KinSrcStep.i"
%include "KinSrcRamp.i"
%include "KinSrcConstRate.i"
%include "KinSrcBrune.i"
%include "KinSrcLiuCos.i"
%include "KinSrcTimeHistory.i"

// End of file
