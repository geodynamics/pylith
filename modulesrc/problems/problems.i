// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

// SWIG interface
%module problems

// Header files for module C++ code
%{
#include "pylith/problems/Problem.hh"
#include "pylith/problems/TimeDependent.hh"
#include "pylith/problems/Physics.hh"
#include "pylith/problems/ObserverSoln.hh"
#include "pylith/problems/ObserverPhysics.hh"
#include "pylith/problems/InitialCondition.hh"
#include "pylith/problems/InitialConditionDomain.hh"
#include "pylith/problems/InitialConditionPatch.hh"
#include "pylith/problems/ProgressMonitor.hh"
#include "pylith/problems/ProgressMonitorTime.hh"
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
%include "../include/physicsarray.i"
%include "../include/outputarray.i"
%include "../include/chararray.i"
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
%include "../utils/PyreComponent.i"
%include "../topology/FieldBase.i"
%include "Problem.i"
%include "TimeDependent.i"
%include "Physics.i"
%include "ObserverSoln.i"
%include "ObserverPhysics.i"
%include "InitialCondition.i"
%include "InitialConditionDomain.i"
%include "InitialConditionPatch.i"
%include "ProgressMonitor.i"
%include "ProgressMonitorTime.i"

// End of file
