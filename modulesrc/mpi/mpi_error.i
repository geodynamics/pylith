// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


// ----------------------------------------------------------------------
// mpi_abort
%inline %{
  void
    mpi_abort(MPI_Comm comm, int err) {
    MPI_Abort(comm, err);
  } // mpi_abort
%}


// End of file
