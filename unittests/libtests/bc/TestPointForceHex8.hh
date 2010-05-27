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
 * @file unittests/libtests/bc/TestPointForceHex8.hh
 *
 * @brief C++ TestPointForce object.
 *
 * C++ unit testing for PointForce for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testpointforcebchex8_hh)
#define pylith_bc_testpointforcebchex8_hh

#include "TestPointForce.hh" // ISA TestPointForce

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForceHex8;
  } // bc
} // pylith

/// C++ unit testing for PointForce for mesh with 2-D tri cells.
class pylith::bc::TestPointForceHex8 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceHex8, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceHex8

#endif // pylith_bc_pointforcebchex8_hh


// End of file 
