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
// destroy_comm
%inline %{
  void
  destroy_comm(MPI_Comm* comm) {
    delete comm;
  } // destroy_comm
%}

// ----------------------------------------------------------------------
// PETSC_COMM_WORLD
%inline %{
  MPI_Comm*
  petsc_comm_world(void) {
    return new MPI_Comm(PETSC_COMM_WORLD);
  } // petc_comm_world
%}

// ----------------------------------------------------------------------
// PETSC_COMM_SELF
%inline %{
  MPI_Comm*
  petsc_comm_self(void) {
    return new MPI_Comm(PETSC_COMM_SELF);
  } // petsc_comm_self
%}

// ----------------------------------------------------------------------
// MPI_COMM_WORLD
%inline %{
  MPI_Comm*
  mpi_comm_world(void) {
    return new MPI_Comm(MPI_COMM_WORLD);
  } // comm_world
%}

// ----------------------------------------------------------------------
// MPI_COMM_SELF
%inline %{
  MPI_Comm*
  mpi_comm_self(void) {
    return new MPI_Comm(MPI_COMM_SELF);
  } // comm_self
%}

// ----------------------------------------------------------------------
// rank()
%inline %{
  int
  rank(MPI_Comm* comm) {
    int value = 0;
    MPI_Comm_rank(*comm, &value);
    return value;
  } // rank
%}

// ----------------------------------------------------------------------
// size()
%inline %{
  int
  size(MPI_Comm* comm) {
    int value = 0;
    MPI_Comm_size(*comm, &value);
    return value;
  } // rank
%}

// ----------------------------------------------------------------------
// barrier()
%inline %{
  void
  barrier(MPI_Comm* comm) {
    MPI_Barrier(*comm);
  } // rank
%}


// End of file

