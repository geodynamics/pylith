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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

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

