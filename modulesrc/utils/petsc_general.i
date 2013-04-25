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

// ----------------------------------------------------------------------
// initialize
%inline %{
  int
  initialize(int argc,
	     char** argv)
  { // initialize
    PetscErrorCode err = 
      PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); CHKERRQ(err);
    return 0;
  } // initialize
%} // inline

// ----------------------------------------------------------------------
// finalize
%inline %{
  int
  finalize(void)
  { // finalize
    PetscErrorCode err = PetscFinalize(); CHKERRQ(err);
    return 0;
  } // finalize
%} // inline

// ----------------------------------------------------------------------
// PetscOptionsSetValue
%inline %{
  int
  optionsSetValue(const char* name,
		  const char* value)
  { // optionsSetValue
    PetscErrorCode err = PetscOptionsSetValue(name, value); CHKERRQ(err);
    return 0;
  } // optionsSetValue
%} // inline

// ----------------------------------------------------------------------
// PetscOptionsHasName
%inline %{
  bool
  optionsHasName(const char* name)
  { // optionsHasName
    PetscBool hasName = PetscBool(0);
    PetscErrorCode err =
      PetscOptionsHasName(PETSC_NULL, name, &hasName);

    return (hasName) ? true : false;
  } // optionsHasName
%} // inline


// End of file

