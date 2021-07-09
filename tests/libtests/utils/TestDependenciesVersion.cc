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

#include "TestDependenciesVersion.hh" // Implementation of class methods

#include "pylith/utils/DependenciesVersion.hh" // USES DependenciesVersion

#include <string.h> // USES strlen()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::utils::TestDependenciesVersion );

// ----------------------------------------------------------------------
// Test mpiVersion()
void
pylith::utils::TestDependenciesVersion::testMPIVersion(void)
{ // testMPIVersion
  const char* v = DependenciesVersion::mpiVersion();
  CPPUNIT_ASSERT(strlen(v) > 0);
} // testMPIVersion

// ----------------------------------------------------------------------
// Test mpiVersion()
void
pylith::utils::TestDependenciesVersion::testMPIImplementation(void)
{ // testMPIImplementation
  const char* v = DependenciesVersion::mpiImplementation();
  CPPUNIT_ASSERT(strlen(v) > 0);
} // testMPIImplementation

// ----------------------------------------------------------------------
// Test mpiStandard()
void
pylith::utils::TestDependenciesVersion::testMPIStandard(void)
{ // testMPIStandard
  const char* v = DependenciesVersion::mpiStandard();
  CPPUNIT_ASSERT(strlen(v) > 0);
} // testMPIStandard

// ----------------------------------------------------------------------
// Test netcdfVersion()
void
pylith::utils::TestDependenciesVersion::testNetCDFVersion(void)
{ // testNetCDFVersion
  const char* v = DependenciesVersion::netcdfVersion();
  CPPUNIT_ASSERT(strlen(v) > 0);
} // testNetCDFVersion

// ----------------------------------------------------------------------
// Test hdf5Version()
void
pylith::utils::TestDependenciesVersion::testHDF5Version(void)
{ // testHDF5Version
  const char* v = DependenciesVersion::hdf5Version();
  CPPUNIT_ASSERT(strlen(v) > 0);
} // testHDF5Version



// End of file 
