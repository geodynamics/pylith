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
 * @file unittests/libtests/bc/TestDirichletPointsMultiTet4.hh
 *
 * @brief C++ TestDirichletPoints object.
 *
 * C++ unit testing for DirichletPoints for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletpointsmultitet4_hh)
#define pylith_bc_testdirichletpointsmultitet4_hh

#include "TestDirichletPointsMulti.hh" // ISA TestDirichletPoints

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletPointsMultiTet4;
  } // bc
} // pylith

/// C++ unit testing for DirichletPoints for mesh with 2-D tri cells.
class pylith::bc::TestDirichletPointsMultiTet4 : public TestDirichletPointsMulti
{ // class TestDirichletPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletPointsMultiTet4 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletPointsMultiTet4

#endif // pylith_bc_dirichletpointsmultitet4_hh


// End of file 
