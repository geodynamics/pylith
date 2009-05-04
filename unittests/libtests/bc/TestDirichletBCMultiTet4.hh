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
 * @file unittests/libtests/bc/TestDirichletBCMultiTet4.hh
 *
 * @brief C++ TestDirichletBC object.
 *
 * C++ unit testing for DirichletBC for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testdirichletbcmultitet4_hh)
#define pylith_bc_testdirichletbcmultitet4_hh

#include "TestDirichletBCMulti.hh" // ISA TestDirichletBC

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestDirichletBCMultiTet4;
  } // bc
} // pylith

/// C++ unit testing for DirichletBC for mesh with 2-D tri cells.
class pylith::bc::TestDirichletBCMultiTet4 : public TestDirichletBCMulti
{ // class TestDirichletBC

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDirichletBCMultiTet4 );
  CPPUNIT_TEST( testSetConstraintSizes );
  CPPUNIT_TEST( testSetConstraints );
  CPPUNIT_TEST( testSetField );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestDirichletBCMultiTet4

#endif // pylith_bc_dirichletbcmultitet4_hh


// End of file 
