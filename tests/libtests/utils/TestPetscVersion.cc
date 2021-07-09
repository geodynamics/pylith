// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPetscVersion.hh" // Implementation of class methods

#include "pylith/utils/PetscVersion.hh" // USES PetscVersion

#include "petsc.h" // USES PETSC_VERSION_*

#include <string> // USES std::string()
#include <string.h> // USES strlen()
#include <stdio.h> // USES snprintf()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestPetscVersion );

// ----------------------------------------------------------------------
// Test isRelease()
void
pylith::utils::TestPetscVersion::testIsRelease(void)
{ // testIsRelease
#if PETSC_VERSION_RELEASE
  CPPUNIT_ASSERT(PetscVersion::isRelease());
#else
  CPPUNIT_ASSERT(!PetscVersion::isRelease());
#endif
} // testIsRelease

// ----------------------------------------------------------------------
// Test version()
void
pylith::utils::TestPetscVersion::testVersion(void)
{ // testVersion
  const int maxsize = 64;
  char value[maxsize];
  snprintf(value, maxsize-1, "%d.%d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR);
  CPPUNIT_ASSERT_EQUAL(std::string(value), std::string(PetscVersion::version()));
} // testVersion

// ----------------------------------------------------------------------
// Test gitRevision()
void
pylith::utils::TestPetscVersion::testGitRevision(void)
{ // testGitRevision
#if PETSC_VERSION_RELEASE
#else
  // Git revision should be of the form vX.X.X-gXXXXX.
  const char* rev = PetscVersion::gitRevision();
  CPPUNIT_ASSERT_EQUAL('v', rev[0]);
#endif
} // testGitRevision

// ----------------------------------------------------------------------
// Test gitDate()
void
pylith::utils::TestPetscVersion::testGitDate(void)
{ // testGitDate
  const char* datetime = PetscVersion::gitDate();
  CPPUNIT_ASSERT(strlen(datetime) > 0);
} // testGitDate

// ----------------------------------------------------------------------
// Test gitBranch()
void
pylith::utils::TestPetscVersion::testGitBranch(void)
{ // testGitBranch
#if PETSC_VERSION_RELEASE
#else
  const char* branch = PetscVersion::gitBranch();
  CPPUNIT_ASSERT(strlen(branch) > 0);
#endif
} // testGitBranch


// ----------------------------------------------------------------------
// Test petscDir()
void
pylith::utils::TestPetscVersion::testPetscDir(void)
{ // testPetscDir
  CPPUNIT_ASSERT(strlen(PetscVersion::petscDir()) > 0);
} // testPetscDir


// ----------------------------------------------------------------------
// Test petscArch()
void
pylith::utils::TestPetscVersion::testPetscArch(void)
{ // testPetscArch
    // Not defined for PETSc prefix installs, so use minimal test.
  CPPUNIT_ASSERT(strlen(PetscVersion::petscArch()) >= 0);
} // testPetscArch


// End of file 
