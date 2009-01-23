// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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

// End of file

