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
%module topology

// Header files for module C++ code
%{
#include "pylith/topology/Mesh.hh"
#include "pylith/topology/MeshOps.hh"
#include "pylith/topology/FieldBase.hh"
#include "pylith/topology/Field.hh"
#include "pylith/topology/Distributor.hh"
#include "pylith/topology/RefineUniform.hh"
#include "pylith/topology/ReverseCuthillMcKee.hh"
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
%include "Mesh.i"
%include "MeshOps.i"
%include "FieldBase.i"
%include "Field.i"
%include "Distributor.i"
%include "RefineUniform.i"
%include "ReverseCuthillMcKee.i"

// End of file

