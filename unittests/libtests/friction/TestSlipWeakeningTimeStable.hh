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
 * @file unittests/libtests/friction/TestSlipWeakeningTimeStable.hh
 *
 * @brief C++ TestSlipWeakeningTimeStable object
 *
 * C++ unit testing for SlipWeakeningTimeStable.
 */

#if !defined(pylith_friction_testslipweakeningtimestable_hh)
#define pylith_friction_testslipweakeningtimestable_hh

#include "TestFrictionModel.hh"

/// Namespace for pylith package
namespace pylith {
  namespace friction {
    class TestSlipWeakeningTimeStable;
  } // friction
} // pylith

/// C++ unit testing for SlipWeakeningTimeStable
class pylith::friction::TestSlipWeakeningTimeStable : public TestFrictionModel
{ // class TestSlipWeakeningTimeStable

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestSlipWeakeningTimeStable );

  CPPUNIT_TEST( testPropertiesMetadata );
  CPPUNIT_TEST( testStateVarsMetadata );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( testHasPropStateVar );
  CPPUNIT_TEST( test_calcFriction );
  CPPUNIT_TEST( test_calcFrictionDeriv );
  CPPUNIT_TEST( test_updateStateVars );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test properties metadata.
  void testPropertiesMetadata(void);

  /// Test state variable metadata.
  void testStateVarsMetadata(void);

  /// Test hasPropStateVar().
  void testHasPropStateVar(void);

}; // class TestSlipWeakeningTimeStable

#endif // pylith_friction_testslipweakeningtimestable_hh


// End of file 
