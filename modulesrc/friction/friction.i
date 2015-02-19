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
%module friction

// Header files for module C++ code
%{
#include "pylith/friction/frictionfwd.hh" // forward declarations

#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "pylith/friction/FrictionModel.hh"
#include "pylith/friction/StaticFriction.hh"
#include "pylith/friction/SlipWeakening.hh"
#include "pylith/friction/SlipWeakeningTime.hh"
#include "pylith/friction/SlipWeakeningTimeStable.hh"
#include "pylith/friction/RateStateAgeing.hh"
#include "pylith/friction/TimeWeakening.hh"

#include "pylith/utils/types.hh"
#include "pylith/utils/array.hh"
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
%include "FrictionModel.i"
%include "StaticFriction.i"
%include "SlipWeakening.i"
%include "SlipWeakeningTime.i"
%include "SlipWeakeningTimeStable.i"
%include "RateStateAgeing.i"
%include "TimeWeakening.i"


// End of file

