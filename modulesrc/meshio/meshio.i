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
%module meshio

// Header files for module C++ code
%{
#include "pylith/meshio/MeshIO.hh"
#include "pylith/meshio/MeshIOAscii.hh"
#include "pylith/meshio/MeshIOLagrit.hh"
#if defined(ENABLE_CUBIT)
#include "pylith/meshio/MeshIOCubit.hh"
#endif

#include "pylith/meshio/FieldFilter.hh"
#include "pylith/meshio/FieldFilterNone.hh"
#include "pylith/meshio/FieldFilterProject.hh"
#include "pylith/meshio/OutputTrigger.hh"
#include "pylith/meshio/OutputTriggerStep.hh"
#include "pylith/meshio/OutputTriggerTime.hh"
#include "pylith/meshio/DataWriter.hh"
#include "pylith/meshio/DataWriterVTK.hh"
#if defined(ENABLE_HDF5)
#include "pylith/meshio/DataWriterHDF5.hh"
#include "pylith/meshio/DataWriterHDF5Ext.hh"
#endif
#include "pylith/meshio/OutputManager.hh"
#include "pylith/meshio/OutputSoln.hh"
  //#include "pylith/meshio/OutputSolnSubset.hh"
  //#include "pylith/meshio/OutputSolnPoints.hh"
#include "pylith/meshio/OutputIntegrator.hh"
#include "pylith/meshio/OutputConstraint.hh"

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

%include "std_string.i"
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
%include "../feassemble/Observer.i"
%include "MeshIOObj.i"
%include "MeshIOAscii.i"
%include "MeshIOLagrit.i"
#if defined(ENABLE_CUBIT)
%include "MeshIOCubit.i"
#endif

%include "FieldFilter.i"
%include "FieldFilterNone.i"
%include "FieldFilterProject.i"
%include "OutputTrigger.i"
%include "OutputTriggerStep.i"
%include "OutputTriggerTime.i"
%include "DataWriter.i"
%include "DataWriterVTK.i"
#if defined(ENABLE_HDF5)
%include "DataWriterHDF5.i"
%include "DataWriterHDF5Ext.i"
#endif
%include "OutputManager.i"
%include "OutputSoln.i"
 //%include "OutputSolnSubset.i"
 //%include "OutputSolnPoints.i"
%include "OutputIntegrator.i"
%include "OutputConstraint.i"


// End of file
