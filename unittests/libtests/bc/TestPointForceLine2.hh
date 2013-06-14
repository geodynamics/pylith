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
 * @file unittests/libtests/bc/TestPointForceLine2.hh
 *
 * @brief C++ TestPointForce object.
 *
 * C++ unit testing for PointForce for mesh with 1-D line cells.
 */

#if !defined(pylith_bc_testpointforcebcline2_hh)
#define pylith_bc_testpointforcebcline2_hh

#include "TestPointForce.hh" // ISA TestPointForce

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestPointForceLine2;
  } // bc
} // pylith

/// C++ unit testing for PointForce for mesh with 2-D tri cells.
class pylith::bc::TestPointForceLine2 : public TestPointForce
{ // class TestPointForce

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUB_SUITE( TestPointForceLine2, TestPointForce );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testVerifyConfiguration );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestPointForceLine2

#endif // pylith_bc_pointforcebcline2_hh


// End of file 
