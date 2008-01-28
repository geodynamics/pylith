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
 * @file unittests/libtests/bc/TestDirichletPointsQuad4.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletpointsquad4_hh)
#define pylith_bc_testdirichletpointsquad4_hh

#include "TestDirichletPoints.hh" // ISA TestDirichletPoints

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsQuad4;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints for mesh with 2-D quad cells.
class pylith::bc::TestDirichletPointsQuad4 : public TestDirichletPoints
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletPointsQuad4, TestDirichletPoints );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletPointsQuad4

#endif // pylith_bc_dirichletpointsquad4_hh


// End of file 
