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

#include "TestDruckerPragerPlaneStrain.hh" // Implementation of class methods

#include "data/DruckerPragerPlaneStrainElasticData.hh" // USES DruckerPragerPlaneStrainElasticData
#include "data/DruckerPragerPlaneStrainTimeDepData.hh" // USES DruckerPragerPlaneStrainTimeDepData

#include "pylith/materials/DruckerPragerPlaneStrain.hh" // USES DruckerPragerPlaneStrain

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestDruckerPragerPlaneStrain );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestDruckerPragerPlaneStrain::setUp(void)
{ // setUp
  _material = new DruckerPragerPlaneStrain();
  DruckerPragerPlaneStrain* matTmp = new DruckerPragerPlaneStrain();
  CPPUNIT_ASSERT(matTmp);
  matTmp->allowTensileYield(true);
  _matElastic = matTmp;

  _data = new DruckerPragerPlaneStrainElasticData();
  _dataElastic = new DruckerPragerPlaneStrainElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestDruckerPragerPlaneStrain::testTimeStep(void)
{ // testTimeStep
  DruckerPragerPlaneStrain material;

  CPPUNIT_ASSERT_EQUAL(false, material._needNewJacobian);

  const PylithScalar dt1 = 1.0;
  material.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(true, material.needNewJacobian());

  const PylithScalar dt2 = 2.0;
  material.timeStep(dt2);
  CPPUNIT_ASSERT_EQUAL(dt2, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(true, material.needNewJacobian());
} // testTimeStep

// ----------------------------------------------------------------------
// Test useElasticBehavior()
void
pylith::materials::TestDruckerPragerPlaneStrain::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  DruckerPragerPlaneStrain material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::DruckerPragerPlaneStrain::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_calcStressElastic ==
      material._calcStressFn);
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastic ==
      material._calcElasticConstsFn);
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastic ==
      material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_calcStressElastoplastic ==
      material._calcStressFn);
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_calcElasticConstsElastoplastic ==
      material._calcElasticConstsFn);
    CPPUNIT_ASSERT(
      &pylith::materials::DruckerPragerPlaneStrain::_updateStateVarsElastoplastic ==
      material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test allowTensileYield()
void
pylith::materials::TestDruckerPragerPlaneStrain::testAllowTensileYield(void)
{ // testAllowTensileYield
  DruckerPragerPlaneStrain material;
  CPPUNIT_ASSERT(!material._allowTensileYield);

  material.allowTensileYield(true);
  CPPUNIT_ASSERT(material._allowTensileYield);
} // testAllowTensileYield

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestDruckerPragerPlaneStrain::testHasStateVars(void)
{ // testHasStateVars
  DruckerPragerPlaneStrain material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPragerPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPragerPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPragerPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestDruckerPragerPlaneStrain::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new DruckerPragerPlaneStrainTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestDruckerPragerPlaneStrain::testHasProperty(void)
{ // testHasProperty
  DruckerPragerPlaneStrain material;

  CPPUNIT_ASSERT(material.hasProperty("mu"));
  CPPUNIT_ASSERT(material.hasProperty("lambda"));
  CPPUNIT_ASSERT(material.hasProperty("density"));
  CPPUNIT_ASSERT(material.hasProperty("alpha_yield"));
  CPPUNIT_ASSERT(material.hasProperty("beta"));
  CPPUNIT_ASSERT(material.hasProperty("alpha_flow"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::materials::TestDruckerPragerPlaneStrain::testHasStateVar(void)
{ // testHasStateVar
  DruckerPragerPlaneStrain material;

  CPPUNIT_ASSERT(material.hasStateVar("stress_zz_initial"));
  CPPUNIT_ASSERT(material.hasStateVar("plastic_strain"));
} // testHasStateVar

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::materials::TestDruckerPragerPlaneStrain::testDBToProperties(void)
{ // testDBToProperties
  CPPUNIT_ASSERT(0 != _material);
  
  TestMaterial::testDBToProperties();

  CPPUNIT_ASSERT_EQUAL(false, _material->isJacobianSymmetric());
} // testDBToProperties


// End of file 
