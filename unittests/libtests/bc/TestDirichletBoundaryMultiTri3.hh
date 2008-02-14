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
 * @file unittests/libtests/bc/TestDirichletBoundaryMultiTri3.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletboundarymultitri3_hh)
#define pylith_bc_testdirichletboundarymultitri3_hh

#include "TestDirichletBoundaryMulti.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryMultiTri3;
  } // bc
} // pylith

/// C++ unit testing for DirichletBoundary for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBoundaryMultiTri3 : public TestDirichletBoundaryMulti
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBoundaryMultiTri3 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryMultiTri3

#endif // pylith_bc_dirichletboundarymultitri3_hh


// End of file 
