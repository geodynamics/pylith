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

#include <portinfo>

#include "TestMaxwellPlaneStrain.hh" // Implementation of class methods

#include "data/MaxwellPlaneStrainElasticData.hh" // USES MaxwellPlaneStrainElasticData
#include "data/MaxwellPlaneStrainTimeDepData.hh" // USES MaxwellPlaneStrainTimeDepData

#include "pylith/materials/MaxwellPlaneStrain.hh" // USES MaxwellPlaneStrain

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellPlaneStrain );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaxwellPlaneStrain::setUp(void)
{ // setUp
  _material = new MaxwellPlaneStrain();
  _matElastic = new MaxwellPlaneStrain();
  _data = new MaxwellPlaneStrainElasticData();
  _dataElastic = new MaxwellPlaneStrainElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestMaxwellPlaneStrain::testTimeStep(void)
{ // testTimeStep
  MaxwellPlaneStrain material;

  CPPUNIT_ASSERT_EQUAL(false, material._needNewJacobian);

  const PylithScalar dt1 = 1.0;
  material.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(false, material.needNewJacobian());

  const PylithScalar dt2 = 2.0;
  material.timeStep(dt2);
  CPPUNIT_ASSERT_EQUAL(dt2, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(true, material.needNewJacobian());
} // testTimeStep

// ----------------------------------------------------------------------
// Test useElasticBehavior()
void
pylith::materials::TestMaxwellPlaneStrain::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  MaxwellPlaneStrain material;

  material.useElasticBehavior(true);
  
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::MaxwellPlaneStrain::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_calcElasticConstsViscoelastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellPlaneStrain::_updateStateVarsViscoelastic ==
		   material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestMaxwellPlaneStrain::testHasStateVars(void)
{ // testHasStateVars
  MaxwellPlaneStrain material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestMaxwellPlaneStrain::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellPlaneStrain::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestMaxwellPlaneStrain::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestMaxwellPlaneStrain::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellPlaneStrain::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestMaxwellPlaneStrain::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
   _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestMaxwellPlaneStrain::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new MaxwellPlaneStrainTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestMaxwellPlaneStrain::testHasProperty(void)
{ // testHasProperty
  MaxwellPlaneStrain material;

  CPPUNIT_ASSERT(material.hasProperty("mu"));
  CPPUNIT_ASSERT(material.hasProperty("lambda"));
  CPPUNIT_ASSERT(material.hasProperty("density"));
  CPPUNIT_ASSERT(material.hasProperty("maxwell_time"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::materials::TestMaxwellPlaneStrain::testHasStateVar(void)
{ // testHasStateVar
  MaxwellPlaneStrain material;

  CPPUNIT_ASSERT(material.hasStateVar("total_strain"));
  CPPUNIT_ASSERT(material.hasStateVar("viscous_strain"));
  CPPUNIT_ASSERT(!material.hasStateVar("stress"));
} // testHasStateVar


// End of file 
