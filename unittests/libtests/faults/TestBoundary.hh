// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestBoundary.hh
 *
 * @brief C++ TestBoundary object
 *
 * C++ unit testing for Fault.
 */

#if !defined(pylith_faults_testboundary_hh)
#define pylith_faults_testboundary_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestBoundary;
    class BoundaryData;
  } // faults
} // pylith

/// C++ unit testing for Fault
class pylith::faults::TestBoundary : public CppUnit::TestFixture
{ // class TestBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestBoundary );

  CPPUNIT_TEST( testCreateBoundaryTri3 );
  CPPUNIT_TEST( testCreateBoundaryTet4 );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test createBoundary() with 2-D triangular element.
  void testCreateBoundaryTri3(void);
  /// Test createBoundary() with 3-D tetrahedral element.
  void testCreateBoundaryTet4(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
public :

  /** Test createBoundary().
   *
   * @param fault Fault for cohesive elements.
   * @param data Cohesive element data.
   */
  void _testCreateBoundary(const BoundaryData& data);

}; // class TestBoundary

#endif // pylith_faults_testboundary_hh

// End of file 
