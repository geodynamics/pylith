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

