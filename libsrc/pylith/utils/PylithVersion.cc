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

#include <portinfo>

#include "PylithVersion.hh" // Implementation of class methods

// ----------------------------------------------------------------------
const bool pylith::utils::PylithVersion::_isRelease = int(PYLITH_RELEASE_VERSION);
const char* pylith::utils::PylithVersion::_version = PYLITH_VERSION;
const char* pylith::utils::PylithVersion::_doi = PYLITH_DOI;
const char* pylith::utils::PylithVersion::_gitBranch = PYLITH_GIT_BRANCH;
const char* pylith::utils::PylithVersion::_gitRevision = PYLITH_GIT_REVISION;
const char* pylith::utils::PylithVersion::_gitDate = PYLITH_GIT_DATE;
const char* pylith::utils::PylithVersion::_gitHash = PYLITH_GIT_HASH;

// ----------------------------------------------------------------------
// Default constructor.
pylith::utils::PylithVersion::PylithVersion(void)
{}

// ----------------------------------------------------------------------
// Default destrictor.
pylith::utils::PylithVersion::~PylithVersion(void)
{}

// ----------------------------------------------------------------------
// Is source from a release?
bool
pylith::utils::PylithVersion::isRelease(void)
{ // isRelease
  return _isRelease;
} // isRelease

// ----------------------------------------------------------------------
// Get version number.
const char*
pylith::utils::PylithVersion::version(void)
{ // version
  return _version;
} // version

// ----------------------------------------------------------------------
// Get DOI.
const char*
pylith::utils::PylithVersion::doi(void)
{ // doi
  return _doi;
} // doi

// ----------------------------------------------------------------------
// Get GIT revision.
const char*
pylith::utils::PylithVersion::gitRevision(void)
{ // gitRevision
  return _gitRevision;
} // gitRevision

// ----------------------------------------------------------------------
// Get GIT hash.
const char*
pylith::utils::PylithVersion::gitHash(void)
{ // gitHash
  return _gitHash;
} // gitHash

// ----------------------------------------------------------------------
// Get date of GIT revision.
const char*
pylith::utils::PylithVersion::gitDate(void)
{ // gitDate
  return _gitDate;
} // gitDate

// ----------------------------------------------------------------------
// Get GIT branch.
const char*
pylith::utils::PylithVersion::gitBranch(void)
{ // gitBranch
  return _gitBranch;
} // gitBranch


// End of file 
