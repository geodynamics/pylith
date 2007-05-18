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
 * @file unittests/libtests/bc/TestDirichletQuad4.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletquad4_hh)
#define pylith_bc_testdirichletquad4_hh

#include "TestDirichlet.hh" // ISA TestDirichlet

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletQuad4;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet for mesh with 2-D quad cells.
class pylith::bc::TestDirichletQuad4 : public TestDirichlet
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletQuad4 );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

}; // class TestDirichletQuad4

#endif // pylith_bc_dirichletquad4_hh


// End of file 
