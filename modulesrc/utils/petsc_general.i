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
  initialize(int argc, char** argv)
  { // initialize
    PetscErrorCode err = 
      PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(err);
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


// End of file

