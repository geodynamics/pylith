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
 * @file unittests/libtests/bc/TestTimeDependentPoints.hh
 *
 * @brief C++ TestTimeDependentPoints object.
 *
 * C++ unit testing for TimeDependentPoints.
 */

#if !defined(pylith_bc_testtimedependentpoints_hh)
#define pylith_bc_testtimedependentpoints_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/bc/bcfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class TestTimeDependentPoints;
  } // bc
} // pylith

/// C++ unit testing for TimeDependentPoints.
class pylith::bc::TestTimeDependentPoints : public CppUnit::TestFixture
{ // class TestTimeDependentPoints

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTimeDependentPoints );

  CPPUNIT_TEST( testBCDOF );
  CPPUNIT_TEST( testGetLabel );
  CPPUNIT_TEST( testQueryDatabases );
  CPPUNIT_TEST( testCalculateValueInitial );
  CPPUNIT_TEST( testCalculateValueRate );
  CPPUNIT_TEST( testCalculateValueChange );
  CPPUNIT_TEST( testCalculateValueChangeTH );
  CPPUNIT_TEST( testCalculateValueAll );
  CPPUNIT_TEST( testCalculateValueIncrInitial );
  CPPUNIT_TEST( testCalculateValueIncrRate );
  CPPUNIT_TEST( testCalculateValueIncrChange );
  CPPUNIT_TEST( testCalculateValueIncrChangeTH );
  CPPUNIT_TEST( testCalculateValueIncrAll );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test bcDOF.
  void testBCDOF(void);

  /// Test _getLabel().
  void testGetLabel(void);

  /// Test _queryDatabases().
  void testQueryDatabases(void);

  /// Test _calculateValue() with initial value.
  void testCalculateValueInitial(void);

  /// Test _calculateValue() with rate.
  void testCalculateValueRate(void);

  /// Test _calculateValue() with temporal change.
  void testCalculateValueChange(void);

  /// Test _calculateValue() with temporal change w/time history.
  void testCalculateValueChangeTH(void);

  /// Test _calculateValue() with initial, rate, and temporal change
  /// w/time history.
  void testCalculateValueAll(void);

  /// Test _calculateValueIncr() with initial value.
  void testCalculateValueIncrInitial(void);

  /// Test _calculateValueIncr() with rate.
  void testCalculateValueIncrRate(void);

  /// Test _calculateValueIncr() with temporal change.
  void testCalculateValueIncrChange(void);

  /// Test _calculateValueIncr() with temporal change w/time history.
  void testCalculateValueIncrChangeTH(void);

  /// Test _calculateValueIncr() with initial, rate, and temporal change
  /// w/time history.
  void testCalculateValueIncrAll(void);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  topology::Mesh* _mesh; ///< Finite-element mesh.
  PointForce* _bc; ///< Point force boundary condition as tester.

}; // class TestTimeDependentPoints

#endif // pylith_bc_timedependentpoints_hh


// End of file 
