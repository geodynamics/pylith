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
 * @file unittests/libtests/bc/TestDirichletHex8.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichlethex8_hh)
#define pylith_bc_testdirichlethex8_hh

#include "TestDirichlet.hh" // ISA TestDirichlet

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletHex8;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet for mesh with 3-D hex cells.
class pylith::bc::TestDirichletHex8 : public TestDirichlet
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletHex8 );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletHex8

#endif // pylith_bc_dirichlethex8_hh


// End of file 
