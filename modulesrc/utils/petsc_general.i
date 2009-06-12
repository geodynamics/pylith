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
// sizeofVoidPtr
%inline %{
  int
  sizeofVoidPtr(void)
  { // sizeofVoidPtr
    return sizeof(void*);
  } // sizeofVoidPtr
%} // inline


// End of file

