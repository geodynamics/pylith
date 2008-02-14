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
 * @file unittests/libtests/bc/TestDirichletBoundaryTet4.hh
 *
 * @brief C++ TestDirichletBoundary object.
 *
 * C++ unit testing for DirichletBoundary for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletboundarytet4_hh)
#define pylith_bc_testdirichletboundaryet4_hh

#include "TestDirichletBoundary.hh" // ISA TestDirichletBoundary

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBoundaryTet4;
  } // bc
} // pylith

/// C++ unit testing for DirichletBoundary for mesh with 3-D tet cells.
class pylith::bc::TestDirichletBoundaryTet4 : public TestDirichletBoundary
{ // class TestDirichletBoundary

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletBoundaryTet4, TestDirichletBoundary );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBoundaryTet4

#endif // pylith_bc_dirichletboundarytet4_hh


// End of file 
