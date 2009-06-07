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
 * @file unittests/libtests/bc/TestPointForceQuad4.hh
 *
 * @brief C++ TestPointForce object.
 *
 * C++ unit testing for PointForce for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testpointforcebcquad4_hh)
#define pylith_bc_testpointforcebcquad4_hh

#include "TestPointForce.hh" // ISA TestPointForce

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForceQuad4;
  } // bc
} // pylith

/// C++ unit testing for PointForce for mesh with 2-D tri cells.
class pylith::bc::TestPointForceQuad4 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceQuad4, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidualAssembled );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceQuad4

#endif // pylith_bc_pointforcebcquad4_hh


// End of file 
