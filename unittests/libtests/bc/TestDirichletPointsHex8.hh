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
 * @file unittests/libtests/bc/TestDirichletPointsHex8.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletpointshex8_hh)
#define pylith_bc_testdirichletpointshex8_hh

#include "TestDirichletPoints.hh" // ISA TestDirichletPoints

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsHex8;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints for mesh with 3-D hex cells.
class pylith::bc::TestDirichletPointsHex8 : public TestDirichletPoints
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletPointsHex8, TestDirichletPoints );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletPointsHex8

#endif // pylith_bc_dirichletpointshex8_hh


// End of file 
