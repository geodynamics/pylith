// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


// SWIG interface
%module mpi

// Header files for module C++ code
%{
#include <petsc.h>
%}

%include "typemaps.i"

// Interfaces
%include "mpi_comm.i"
%include "mpi_error.i"
%include "mpi_reduce.i"


// End of file
