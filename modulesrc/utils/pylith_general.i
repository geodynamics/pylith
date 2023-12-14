// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

// ----------------------------------------------------------------------
// sizeofVoidPtr
%inline %{
  int
  sizeofVoidPtr(void)
  { // sizeofVoidPtr
    return sizeof(void*);
  } // sizeofVoidPtr
%} // inline

// ----------------------------------------------------------------------
// sizeofPylithScalar
%inline %{
  int
  sizeofPylithScalar(void)
  { // sizeofPylithScalar
    return sizeof(PylithScalar);
  } // sizeofPylithScalar
%} // inline


// ----------------------------------------------------------------------
// isCUDAEnabled
%inline %{
  bool
  isCUDAEnabled(void)
  { // isCUDAEnabled
#if ENABLE_CUDA
    return true;
#else
    return false;
#endif
  } // isCUDAEnabled
%} // inline


// End of file

