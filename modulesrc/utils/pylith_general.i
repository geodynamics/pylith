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

