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
 * @file unittests/libtests/bc/TestDirichletBoundaryMultiTet4.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletboundarymultitet4_hh)
#define pylith_bc_testdirichletboundarymultitet4_hh

#include "TestDirichletBoundaryMulti.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryMultiTet4;
  } // bc
} // pylith

/// C++ unit testing for DirichletBoundary for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBoundaryMultiTet4 : public TestDirichletBoundaryMulti
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBoundaryMultiTet4 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryMultiTet4

#endif // pylith_bc_dirichletboundarymultitet4_hh


// End of file 
