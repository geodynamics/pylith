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

#include "PetscVersion.hh" // Implementation of class methods

#include "petsc.h"

// ----------------------------------------------------------------------
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define PYLITH_PETSC_VERSION STR(PETSC_VERSION_MAJOR) "." STR(PETSC_VERSION_MINOR) "." STR(PETSC_VERSION_SUBMINOR)
const bool pylith::utils::PetscVersion::_isRelease = PETSC_VERSION_RELEASE;
const char* pylith::utils::PetscVersion::_version = PYLITH_PETSC_VERSION;
#if defined(PETSC_VERSION_BRANCH_GIT)
const char* pylith::utils::PetscVersion::_gitBranch = PETSC_VERSION_BRANCH_GIT;
#else
const char* pylith::utils::PetscVersion::_gitBranch = "branch-not-available";
#endif
const char* pylith::utils::PetscVersion::_gitRevision = PETSC_VERSION_GIT;
const char* pylith::utils::PetscVersion::_gitDate = PETSC_VERSION_DATE_GIT;
const char* pylith::utils::PetscVersion::_petscDir = PETSC_DIR;
const char* pylith::utils::PetscVersion::_petscArch = PETSC_ARCH;

// ----------------------------------------------------------------------
// Default constructor.
pylith::utils::PetscVersion::PetscVersion(void)
{}

// ----------------------------------------------------------------------
// Default destrictor.
pylith::utils::PetscVersion::~PetscVersion(void)
{}

// ----------------------------------------------------------------------
// Is source from a release?
bool
pylith::utils::PetscVersion::isRelease(void)
{ // isRelease
  return _isRelease;
} // isRelease

// ----------------------------------------------------------------------
// Get version number.
const char*
pylith::utils::PetscVersion::version(void)
{ // version
  return _version;
} // version

// ----------------------------------------------------------------------
// Get GIT revision.
const char*
pylith::utils::PetscVersion::gitRevision(void)
{ // gitRevision
  return _gitRevision;
} // gitRevision

// ----------------------------------------------------------------------
// Get date of GIT revision.
const char*
pylith::utils::PetscVersion::gitDate(void)
{ // gitDate
  return _gitDate;
} // gitDate

// ----------------------------------------------------------------------
// Get GIT branch.
const char*
pylith::utils::PetscVersion::gitBranch(void)
{ // gitBranch
  return _gitBranch;
} // gitBranch

// ----------------------------------------------------------------------
// Get PETSC_DIR.
const char*
pylith::utils::PetscVersion::petscDir(void)
{ // petscDir
  return _petscDir;
} // petscDir
  
// ----------------------------------------------------------------------
// Get PETSC_ARCH.
const char*
pylith::utils::PetscVersion::petscArch(void)
{ // petscArch
  return _petscArch;
} // petscArch
  
// End of file 
