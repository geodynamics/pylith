// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/utils/constdefs.i
 *
 * @brief PyLith constants.
 */

// ----------------------------------------------------------------------
// PYLITH_MAXDOUBLE
%inline %{
  double
  maxdouble(void)
  { // maxdouble
    return pylith::PYLITH_MAXDOUBLE;
  } // maxdouble
%} // inline


// ----------------------------------------------------------------------
// PYLITH_MAXFLOAT
%inline %{
  float
  maxfloat(void)
  { // maxfloat
    return pylith::PYLITH_MAXFLOAT;
  } // maxfloat
%} // inline


// ----------------------------------------------------------------------
// PYLITH_MAXSCALAR
%inline %{
  double
  maxscalar(void)
  { // maxscalar
    return pylith::PYLITH_MAXSCALAR;
  } // maxscalar
%} // inline


// End of file 
