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
 * @file unittests/libtests/bc/TestNeumannTri3.hh
 *
 * @brief C++ TestNeumann object.
 *
 * C++ unit testing for Neumann for mesh with 2-D tri cells.
 */

#if !defined(pylith_bc_testneumanntri3_hh)
#define pylith_bc_testneumanntri3_hh

#include "TestNeumann.hh" // ISA TestNeumann

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestNeumannTri3;
  } // bc
} // pylith

/// C++ unit testing for Neumann for mesh with 2-D tri cells.
class pylith::bc::TestNeumannTri3 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestNeumannTri3, TestNeumann );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannTri3

#endif // pylith_bc_neumanntri3_hh


// End of file 
