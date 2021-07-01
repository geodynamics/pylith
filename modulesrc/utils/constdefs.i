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
