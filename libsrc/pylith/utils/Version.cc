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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Version.hh" // Implementation of class methods

// ----------------------------------------------------------------------
// Default constructor.
pylith::utils::Version::Version(void)
{}

// ----------------------------------------------------------------------
// Default destrictor.
pylith::utils::Version::~Version(void)
{}

// ----------------------------------------------------------------------
// Get PyLith releave version number.
const char*
pylith::utils::Version::version(void)
{ // version
  const char* value = PACKAGE_VERSION;

  return value;
} // version


// End of file 
