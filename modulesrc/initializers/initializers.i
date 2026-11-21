// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


// SWIG interface
%module initializers

// Header files for module C++ code
%{
#include "pylith/initializers/Initializer.hh"
#include "pylith/initializers/MeshReader.hh"
#include "pylith/initializers/MeshWriter.hh"
#include "pylith/initializers/MeshReordering.hh"
#include "pylith/initializers/MeshRefiner.hh"
#include "pylith/initializers/MeshDistributor.hh"
#include "pylith/initializers/MeshInsertInterfaces.hh"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
} // exception

%include "phasesarray.i"

%include "Initializer.i"
%include "InitializePhase.i"
%include "MeshReader.i"
%include "MeshWriter.i"
%include "MeshReordering.i"
%include "MeshRefiner.i"
%include "MeshDistributor.i"
%include "MeshInsertInterfaces.i"


// End of file
