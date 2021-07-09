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

/** @file tests/libtests/utils/TestDependenciesVersion.hh
 *
 * @brief C++ TestDependenciesVersion object
 *
 * C++ unit testing for DependenciesVersion.
 */

#if !defined(pylith_testdependenciesversion_hh)
#define pylith_testdependenciesversion_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace utils {
    class TestDependenciesVersion;
  } // utils
} // pylith

/// C++ unit testing for DependenciesVersion
class pylith::utils::TestDependenciesVersion : public CppUnit::TestFixture
{ // class TestDependenciesVersion

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDependenciesVersion );

  CPPUNIT_TEST( testMPIVersion );
  CPPUNIT_TEST( testMPIImplementation );
  CPPUNIT_TEST( testMPIStandard );
  CPPUNIT_TEST( testNetCDFVersion );
  CPPUNIT_TEST( testHDF5Version );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test mpiVersion()
  void testMPIVersion(void);

  /// Test mpiImplementation()
  void testMPIImplementation(void);

  /// Test mpiStandard()
  void testMPIStandard(void);

  /// Test netcdfVersion()
  void testNetCDFVersion(void);

  /// Test hdf5Version()
  void testHDF5Version(void);

}; // class TestDependenciesVersion

#endif // pylith_utils_testdependenciesversion_hh

// End of file 
