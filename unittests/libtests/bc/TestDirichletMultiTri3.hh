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
 * @file unittests/libtests/bc/TestDirichletMultiTri3.hh
 *
 * @brief C++ TestDirichlet object.
 *
 * C++ unit testing for Dirichlet for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletmultitri3_hh)
#define pylith_bc_testdirichletmultitri3_hh

#include "TestDirichletMulti.hh" // ISA TestDirichlet

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletMultiTri3;
  } // bc
} // pylith

/// C++ unit testing for Dirichlet for mesh with 2-D tri cells.
class pylith::bc::TestDirichletMultiTri3 : public TestDirichletMulti
{ // class TestDirichlet

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletMultiTri3 );
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

}; // class TestDirichletMultiTri3

#endif // pylith_bc_dirichletmultitri3_hh


// End of file 
