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
 * @file unittests/libtests/bc/TestNeumannLine2.hh
 *
 * @brief C++ TestNeumann object.
 *
 * C++ unit testing for Neumann for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testneumannline2_hh)
#define pylith_bc_testneumannline2_hh

#include "TestNeumann.hh" // ISA TestNeumann

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestNeumannLine2;
  } // bc
} // pylith

/// C++ unit testing for Neumann for mesh with 1-D line cells.
class pylith::bc::TestNeumannLine2 : public TestNeumann
{ // class TestNeumann

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestNeumannLine2, TestNeumann );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestNeumannLine2

#endif // pylith_bc_neumannline2_hh


// End of file 
