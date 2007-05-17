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
 * @file unittests/libtests/bc/TestBoundaryCondition.hh
 *
 * @brief C++ TestBoundaryCondition object
 *
 * C++ unit testing for BoundaryCondition.
 */

#if !defined(pylith_bc_testboundarycondition_hh)
#define pylith_bc_testboundarycondition_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryCondition;
  } // bc
} // pylith

/// C++ unit testing for BoundaryCondition
class pylith::bc::TestBoundaryCondition : public CppUnit::TestFixture
{ // class TestBoundaryCondition

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestBoundaryCondition );
  CPPUNIT_TEST( testID );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testDB );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test id()
  void testID(void);

  /// Test label()
  void testLabel(void);

  /// Test db()
  void testDB(void);

}; // class TestBoundaryCondition

#endif // pylith_bc_testboundarycondition_hh

// End of file 
