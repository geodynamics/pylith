// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
// SWIG interface
%module problems

// Header files for module C++ code
%{
#include "pylith/problems/Problem.hh"
#include "pylith/problems/TimeDependent.hh"
#include "pylith/problems/GreensFns.hh"
#include "pylith/problems/Physics.hh"
#include "pylith/problems/Observer.hh"
#include "pylith/problems/ObserverSoln.hh"
#include "pylith/problems/ObserverPhysics.hh"
#include "pylith/problems/InitialCondition.hh"
#include "pylith/problems/InitialConditionDomain.hh"
#include "pylith/problems/InitialConditionPatch.hh"
#include "pylith/problems/ProgressMonitor.hh"
#include "pylith/problems/ProgressMonitorTime.hh"
#include "pylith/problems/ProgressMonitorStep.hh"
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
%include "GreensFns.i"
%include "Physics.i"
%include "Observer.i"
%include "ObserverSoln.i"
%include "ObserverPhysics.i"
%include "InitialCondition.i"
%include "InitialConditionDomain.i"
%include "InitialConditionPatch.i"
%include "ProgressMonitor.i"
%include "ProgressMonitorTime.i"
%include "ProgressMonitorStep.i"

// End of file
