// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

// ----------------------------------------------------------------------
// mpi_abort
%inline %{
  void
    mpi_abort(MPI_Comm comm, int err) {
    MPI_Abort(comm, err);
  } // mpi_abort
%}


// End of file

