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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/materials/TestDruckerPragerPlaneStrain.hh
 *
 * @brief C++ TestDruckerPragerPlaneStrain object
 *
 * C++ unit testing for DruckerPragerPlaneStrain.
 */

#if !defined(pylith_materials_testdruckerpragerplanestrain_hh)
#define pylith_materials_testdruckerpragerplanestrain_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestDruckerPragerPlaneStrain;
    class DruckerPragerPlaneStrainElasticData;
    class DruckerPragerPlaneStrainTimeDepData;
  } // materials
} // pylith

/// C++ unit testing for DruckerPragerPlaneStrain
class pylith::materials::TestDruckerPragerPlaneStrain : public TestElasticMaterial
{ // class TestDruckerPragerPlaneStrain

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestDruckerPragerPlaneStrain );

  CPPUNIT_TEST( testDimension );
  CPPUNIT_TEST( testTensorSize );
  CPPUNIT_TEST( testDBToProperties );
  CPPUNIT_TEST( testNonDimProperties );
  CPPUNIT_TEST( testDimProperties );
  CPPUNIT_TEST( testDBToStateVars );
  CPPUNIT_TEST( testNonDimStateVars );
  CPPUNIT_TEST( testDimStateVars );
  CPPUNIT_TEST( test_calcDensity );
  CPPUNIT_TEST( test_stableTimeStepImplicit );

  // Need to test Drucker-Prager elastoplastic specific behavior.
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testUseLinearBehavior );
  CPPUNIT_TEST( testAllowTensileYield );
  CPPUNIT_TEST( testHasStateVars );

  CPPUNIT_TEST( test_calcStressElastic );
  CPPUNIT_TEST( test_calcStressTimeDep );
  CPPUNIT_TEST( test_calcElasticConstsElastic );
  CPPUNIT_TEST( test_calcElasticConstsTimeDep );
  CPPUNIT_TEST( test_updateStateVarsElastic );
  CPPUNIT_TEST( test_updateStateVarsTimeDep );

  CPPUNIT_TEST( testHasProperty );
  CPPUNIT_TEST( testHasStateVar );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

  /// Test timeStep()
  void testTimeStep(void);

  /// Test useLinearBehavior()
  void testUseLinearBehavior(void);

  /// Test allowTensileYield()
  void testAllowTensileYield(void);

  /// Test hasStateVars()
  void testHasStateVars(void);

  /// Test _calcStressElastic()
  void test_calcStressElastic(void);

  /// Test _calcElasticConstsElastic()
  void test_calcElasticConstsElastic(void);

  /// Test _updateStateVarsElastic()
  void test_updateStateVarsElastic(void);

  /// Test _calcStressTimeDep()
  void test_calcStressTimeDep(void);

  /// Test _calcElasticConstsTimeDep()
  void test_calcElasticConstsTimeDep(void);

  /// Test _updateStatevarsTimeDep()
  void test_updateStateVarsTimeDep(void);

  /// Test _stableTimeStepImplicit()
  void test_stableTimeStepImplicit(void);

  /// Test hasProperty()
  void testHasProperty(void);

  /// Test hasStateVar()
  void testHasStateVar(void);

  /// Test _dbToProperties() and check flag for symmetry of Jacobian.
  void testDBToProperties(void);

}; // class TestDruckerPragerPlaneStrain

#endif // pylith_materials_testdruckerpragerplanestrain_hh


// End of file 
