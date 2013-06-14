// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestPointForceTri3.hh
 *
 * @brief C++ TestPointForce object.
 *
 * C++ unit testing for PointForce for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testpointforcebctri3_hh)
#define pylith_bc_testpointforcebctri3_hh

#include "TestPointForce.hh" // ISA TestPointForce

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForceTri3;
  } // bc
} // pylith

/// C++ unit testing for PointForce for mesh with 2-D tri cells.
class pylith::bc::TestPointForceTri3 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceTri3, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceTri3

#endif // pylith_bc_pointforcebctri3_hh


// End of file 
