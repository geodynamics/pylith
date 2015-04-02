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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/bc/TestBoundaryConditionPoints.hh
 *
 * @brief C++ TestBoundaryConditionPoints object.
 *
 * C++ unit testing for BoundaryConditionPoints.
 */

#if !defined(pylith_bc_testboundaryconditionpoints_hh)
#define pylith_bc_testboundaryconditionpoints_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestBoundaryConditionPoints;
    class PointForceData;
  } // bc
} // pylith

/// C++ unit testing for PointForce.
class pylith::bc::TestBoundaryConditionPoints : public CppUnit::TestFixture
{ // class TestBoundaryConditionPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestBoundaryConditionPoints );

  CPPUNIT_TEST( testGetPoints );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test _getPoints().
  void testGetPoints(void);

}; // class TestBoundaryConditionPoints

#endif // pylith_bc_testboundaryconditionpoints_hh


// End of file 
