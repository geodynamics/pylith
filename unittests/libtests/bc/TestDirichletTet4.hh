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
 * @file unittests/libtests/bc/TestDirichletTet4.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichlettet4_hh)
#define pylith_bc_testdirichletet4_hh

#include "TestDirichlet.hh" // ISA TestDirichlet

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletTet4;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet for mesh with 3-D tet cells.
class pylith::bc::TestDirichletTet4 : public TestDirichlet
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletTet4 );
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

}; // class TestDirichletTet4

#endif // pylith_bc_dirichlettet4_hh


// End of file 
