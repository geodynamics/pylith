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
 * @file unittests/libtests/materials/TestPowerLaw3D.hh
 *
 * @brief C++ TestPowerLaw3D object
 *
 * C++ unit testing for PowerLaw3D.
 */

#if !defined(pylith_materials_testpowerlaw3d_hh)
#define pylith_materials_testpowerlaw3d_hh

#include "TestElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class TestPowerLaw3D;
    class PowerLaw3DElasticData;
    class PowerLaw3DTimeDepData;
  } // materials
} // pylith

/// C++ unit testing for PowerLaw3D
class pylith::materials::TestPowerLaw3D : public TestElasticMaterial
{ // class TestPowerLaw3D

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestPowerLaw3D );

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
  CPPUNIT_TEST( test_stableTimeStepExplicit );

  // Need to test Power Law viscoelastic specific behavior.
  CPPUNIT_TEST( testTimeStep );
  CPPUNIT_TEST( testUseElasticBehavior );
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

  /// Test useElasticBehavior()
  void testUseElasticBehavior(void);

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

}; // class TestPowerLaw3D

#endif // pylith_materials_testpowerlaw3d_hh


// End of file 
