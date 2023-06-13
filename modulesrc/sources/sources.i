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
%module sources

// Header files for module C++ code
%{
#include "pylith/sources/Source.hh"
#include "pylith/sources/WellboreSource.hh"
#include "pylith/sources/SquarePulseSource.hh"
#include "pylith/sources/PointForce.hh"
#include "pylith/sources/MomentTensorForce.hh"
#include "pylith/sources/RickerWavelet.hh"
#include "pylith/sources/GaussianWavelet.hh"
#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh"

#include "pylith/utils/arrayfwd.hh"
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

%include "Source.i"
%include "WellboreSource.i"
%include "SquarePulseSource.i"
%include "PointForce.i"
%include "MomentTensorForce.i"
%include "SourceTimeFunctionMomentTensorForce.i"
%include "RickerWavelet.i"
%include "GaussianWavelet.i"

// End of file
