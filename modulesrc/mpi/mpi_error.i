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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
// mpi_abort
%inline %{
  void
    mpi_abort(MPI_Comm comm, int err) {
    MPI_Abort(comm, err);
  } // mpi_abort
%}


// End of file

