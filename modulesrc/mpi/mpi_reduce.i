// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

// ----------------------------------------------------------------------
// MPI_SUM
%inline %{
  MPI_Op*
  mpi_sum(void) {
    return new MPI_Op(MPI_SUM);
  } // mpi_sum
%}


// ----------------------------------------------------------------------
// MPI_MIN
%inline %{
  MPI_Op*
  mpi_min(void) {
    return new MPI_Op(MPI_MIN);
  } // mpi_min
%}


// ----------------------------------------------------------------------
// MPI_MAX
%inline %{
  MPI_Op*
  mpi_max(void) {
    return new MPI_Op(MPI_MAX);
  } // mpi_max
%}


// ----------------------------------------------------------------------
// allreduce_scalar_double
%inline %{
  double
    allreduce_scalar_double(double value,
			    MPI_Op* op,
			    MPI_Comm* comm) {
    double result = 0.0;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, *op, *comm);
    return result;
  } // allreduce_scalar_double
%}


// ----------------------------------------------------------------------
// allreduce_scalar_int
%inline %{
  int
    allreduce_scalar_int(int value,
			 MPI_Op* op,
			 MPI_Comm* comm) {
    int result = 0;
    MPI_Allreduce(&value, &result, 1, MPI_INT, *op, *comm);
    return result;
  } // allreduce_int
%}


// End of file

