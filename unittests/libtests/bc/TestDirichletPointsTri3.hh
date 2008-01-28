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
 * @file unittests/libtests/bc/TestDirichletPointsTri3.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletpointstri3_hh)
#define pylith_bc_testdirichletpointstri3_hh

#include "TestDirichletPoints.hh" // ISA TestDirichletPoints

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsTri3;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints for mesh with 2-D tri cells.
class pylith::bc::TestDirichletPointsTri3 : public TestDirichletPoints
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestDirichletPointsTri3, TestDirichletPoints );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletPointsTri3

#endif // pylith_bc_dirichletpointstri3_hh


// End of file 
