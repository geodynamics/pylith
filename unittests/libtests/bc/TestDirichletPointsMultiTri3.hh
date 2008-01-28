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
 * @file unittests/libtests/bc/TestDirichletPointsMultiTri3.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletpointsmultitri3_hh)
#define pylith_bc_testdirichletpointsmultitri3_hh

#include "TestDirichletPointsMulti.hh" // ISA TestDirichletPoints

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsMultiTri3;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints for mesh with 2-D tri cells.
class pylith::bc::TestDirichletPointsMultiTri3 : public TestDirichletPointsMulti
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletPointsMultiTri3 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletPointsMultiTri3

#endif // pylith_bc_dirichletpointsmultitri3_hh


// End of file 
