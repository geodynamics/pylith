// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/PetscVersion.hh" // Implementation of class methods

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
pylith::utils::PetscVersion::PetscVersion(void) {}


// ----------------------------------------------------------------------
// Default destrictor.
pylith::utils::PetscVersion::~PetscVersion(void) {}


// ----------------------------------------------------------------------
// Is source from a release?
bool
pylith::utils::PetscVersion::isRelease(void) { // isRelease
    return _isRelease;
} // isRelease


// ----------------------------------------------------------------------
// Get version number.
const char*
pylith::utils::PetscVersion::version(void) { // version
    return _version;
} // version


// ----------------------------------------------------------------------
// Get GIT revision.
const char*
pylith::utils::PetscVersion::gitRevision(void) { // gitRevision
    return _gitRevision;
} // gitRevision


// ----------------------------------------------------------------------
// Get date of GIT revision.
const char*
pylith::utils::PetscVersion::gitDate(void) { // gitDate
    return _gitDate;
} // gitDate


// ----------------------------------------------------------------------
// Get GIT branch.
const char*
pylith::utils::PetscVersion::gitBranch(void) { // gitBranch
    return _gitBranch;
} // gitBranch


// ----------------------------------------------------------------------
// Get PETSC_DIR.
const char*
pylith::utils::PetscVersion::petscDir(void) { // petscDir
    return _petscDir;
} // petscDir


// ----------------------------------------------------------------------
// Get PETSC_ARCH.
const char*
pylith::utils::PetscVersion::petscArch(void) { // petscArch
    return _petscArch;
} // petscArch


// End of file
