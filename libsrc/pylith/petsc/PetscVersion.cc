// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "PetscVersion.hh" // Implementation of class methods

#include "petscsys.h"

// ----------------------------------------------------------------------
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define PYLITH_PETSC_VERSION STR(PETSC_VERSION_MAJOR) "." STR(PETSC_VERSION_MINOR) "." STR(PETSC_VERSION_SUBMINOR)
const bool pylith::petsc::PetscVersion::_isRelease = PETSC_VERSION_RELEASE;
const char* pylith::petsc::PetscVersion::_version = PYLITH_PETSC_VERSION;
#if defined(PETSC_VERSION_BRANCH_GIT)
const char* pylith::petsc::PetscVersion::_gitBranch = PETSC_VERSION_BRANCH_GIT;
#else
const char* pylith::petsc::PetscVersion::_gitBranch = "branch-not-available";
#endif
const char* pylith::petsc::PetscVersion::_gitRevision = PETSC_VERSION_GIT;
const char* pylith::petsc::PetscVersion::_gitDate = PETSC_VERSION_DATE_GIT;
const char* pylith::petsc::PetscVersion::_petscDir = PETSC_DIR;
const char* pylith::petsc::PetscVersion::_petscArch = PETSC_ARCH;

// ----------------------------------------------------------------------
// Default constructor.
pylith::petsc::PetscVersion::PetscVersion(void) {}


// ----------------------------------------------------------------------
// Default destructor.
pylith::petsc::PetscVersion::~PetscVersion(void) {}


// ----------------------------------------------------------------------
// Is source from a release?
bool
pylith::petsc::PetscVersion::isRelease(void) { // isRelease
    return _isRelease;
} // isRelease


// ----------------------------------------------------------------------
// Get version number.
const char*
pylith::petsc::PetscVersion::version(void) { // version
    return _version;
} // version


// ----------------------------------------------------------------------
// Get GIT revision.
const char*
pylith::petsc::PetscVersion::gitRevision(void) { // gitRevision
    return _gitRevision;
} // gitRevision


// ----------------------------------------------------------------------
// Get date of GIT revision.
const char*
pylith::petsc::PetscVersion::gitDate(void) { // gitDate
    return _gitDate;
} // gitDate


// ----------------------------------------------------------------------
// Get GIT branch.
const char*
pylith::petsc::PetscVersion::gitBranch(void) { // gitBranch
    return _gitBranch;
} // gitBranch


// ----------------------------------------------------------------------
// Get PETSC_DIR.
const char*
pylith::petsc::PetscVersion::petscDir(void) { // petscDir
    return _petscDir;
} // petscDir


// ----------------------------------------------------------------------
// Get PETSC_ARCH.
const char*
pylith::petsc::PetscVersion::petscArch(void) { // petscArch
    return _petscArch;
} // petscArch


// End of file
