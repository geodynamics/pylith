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
// initialize
%inline %{
  int
  initialize(int argc,
	     char** argv)
  { // initialize
    PylithCallPetsc(PetscInitialize(&argc, &argv, NULL, NULL));
    return 0;
  } // initialize
%} // inline

// ----------------------------------------------------------------------
// finalize
%inline %{
  int
  finalize(void)
  { // finalize
    PylithCallPetsc(PetscFinalize());
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
    PylithCallPetsc(PetscOptionsSetValue(NULL, name, value));
    return 0;
  } // optionsSetValue
%} // inline

// ----------------------------------------------------------------------
// PetscOptionsHasName
%inline %{
  bool
    optionsHasName(const char* name,
		   const char* pre)
  { // optionsHasName
    PetscBool hasName = PetscBool(0);
    PylithCallPetsc(PetscOptionsHasName(NULL, pre, name, &hasName));

    return (hasName) ? true : false;
  } // optionsHasName
%} // inline


// ----------------------------------------------------------------------
// PetscCitationsRegister
%inline %{
  int
    citationsRegister(const char* entry)
  { // citationsRegister
    PetscBool set = PetscBool(PETSC_FALSE);
    PylithCallPetsc(PetscCitationsRegister(entry, &set));

    return 0;
  } // citationsRegister
%} // inline


// End of file

