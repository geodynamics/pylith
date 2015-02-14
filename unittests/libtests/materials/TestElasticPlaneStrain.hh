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
 * @file unittests/libtests/materials/TestElasticPlaneStrain.hh
 *
 * @brief C++ TestElasticPlaneStrain object
 *
 * C++ unit testing for ElasticPlaneStrain.
 */

#if !defined(pylith_materials_testelasticplanestrain_hh)
#define pylith_materials_testelasticplanestrain_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestElasticPlaneStrain;
  } // materials
} // pylith

/// C++ unit testing for ElasticPlaneStrain
class pylith::materials::TestElasticPlaneStrain : public TestElasticMaterial
{ // class TestElasticPlaneStrain

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestElasticPlaneStrain );

  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testTensorSize );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( testStableTimeStepImplicit );
  CPPUNIT_TEST( test_calcDensity );
  CPPUNIT_TEST( test_calcStress );
  CPPUNIT_TEST( test_calcElasticConsts );
  CPPUNIT_TEST( test_updateStateVars );
  CPPUNIT_TEST( test_stableTimeStepImplicit );
  CPPUNIT_TEST( test_stableTimeStepExplicit );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test stableTimeStepImplicit().
  void testStableTimeStepImplicit(void);

}; // class TestElasticPlaneStrain

#endif // pylith_materials_testelasticplanestrain_hh


// End of file 
