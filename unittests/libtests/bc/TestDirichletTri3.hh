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
 * @file unittests/libtests/bc/TestDirichletTri3.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichlettri3_hh)
#define pylith_bc_testdirichlettri3_hh

#include "TestDirichlet.hh" // ISA TestDirichlet

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletTri3;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet for mesh with 2-D tri cells.
class pylith::bc::TestDirichletTri3 : public TestDirichlet
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletTri3 );
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

}; // class TestDirichletTri3

#endif // pylith_bc_dirichlettri3_hh


// End of file 
