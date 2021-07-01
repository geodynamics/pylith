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
// mat_assemble
// Assemble matrix.
%inline %{
  int
  mat_assemble(Mat* mat,
               const char* mode)
  { // mat_assemble
    PetscErrorCode err = 0;
    if (0 == strcmp(mode, "final_assembly")) {
      err = MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);CHKERRQ(err);
      err = MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);CHKERRQ(err);
    } else if (0 == strcmp(mode, "flush_assembly")) {
      err = MatAssemblyBegin(*mat, MAT_FLUSH_ASSEMBLY); CHKERRQ(err);
      err = MatAssemblyEnd(*mat, MAT_FLUSH_ASSEMBLY); CHKERRQ(err);
    } else
      throw std::runtime_error("Unknown mode");
  return 0;
  } // mat_assemble
%} // inline

// ----------------------------------------------------------------------
// mat_setzero
// Zero out entries in matrix (retain structure).
%inline %{
  int
  mat_setzero(Mat* mat)
  { // mat_setzero
    PetscErrorCode err = MatZeroEntries(*mat); CHKERRQ(err);
    return 0;
  } // mat_setzero
%} // inline

// ----------------------------------------------------------------------
// mat_view
// View matrix.
%inline %{
  int
  mat_view(Mat* mat)
  { // mat_view
    PetscErrorCode err = 
      MatView(*mat, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(err);
    return 0;
  } // mat_view
%} // inline

// ----------------------------------------------------------------------
// mat_view_binary
// Write matrix to binary file.
%inline %{
  int
  mat_view_binary(Mat* mat,
		    const char* filename)
  { // mat_view_binary
  PetscViewer viewer;
  PetscErrorCode err = 
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename,
			  FILE_MODE_WRITE, &viewer); CHKERRQ(err);
  err = MatView(*mat, viewer); CHKERRQ(err);
  err = PetscViewerDestroy(viewer); CHKERRQ(err);
  return 0;
  } // mat_view_binary
%} // inline


// End of file

