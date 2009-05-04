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
 * @file unittests/libtests/bc/TestDirichletBCMultiTri3.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * C++ unit testing for DirichletBC for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletbcmultitri3_hh)
#define pylith_bc_testdirichletbcmultitri3_hh

#include "TestDirichletBCMulti.hh" // ISA TestDirichletBC

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCMultiTri3;
  } // bc
} // pylith

/// C++ unit testing for DirichletBC for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBCMultiTri3 : public TestDirichletBCMulti
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBCMultiTri3 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCMultiTri3

#endif // pylith_bc_dirichletbcmultitri3_hh


// End of file 
