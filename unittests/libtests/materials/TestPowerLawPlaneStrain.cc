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

#include "TestPowerLawPlaneStrain.hh" // Implementation of class methods

#include "data/PowerLawPlaneStrainElasticData.hh" // USES PowerLawPlaneStrainElasticData
#include "data/PowerLawPlaneStrainTimeDepData.hh" // USES PowerLawPlaneStrainTimeDepData

#include "pylith/materials/PowerLawPlaneStrain.hh" // USES PowerLawPlaneStrain

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestPowerLawPlaneStrain );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestPowerLawPlaneStrain::setUp(void)
{ // setUp
  _material = new PowerLawPlaneStrain();
  _matElastic = new PowerLawPlaneStrain();
  _data = new PowerLawPlaneStrainElasticData();
  _dataElastic = new PowerLawPlaneStrainElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestPowerLawPlaneStrain::testTimeStep(void)
{ // testTimeStep
  PowerLawPlaneStrain material;

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
pylith::materials::TestPowerLawPlaneStrain::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  PowerLawPlaneStrain material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::PowerLawPlaneStrain::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_calcElasticConstsViscoelastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLawPlaneStrain::_updateStateVarsViscoelastic ==
		   material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestPowerLawPlaneStrain::testHasStateVars(void)
{ // testHasStateVars
  PowerLawPlaneStrain material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestPowerLawPlaneStrain::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestPowerLawPlaneStrain::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestPowerLawPlaneStrain::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestPowerLawPlaneStrain::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLawPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestPowerLawPlaneStrain::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLawPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestPowerLawPlaneStrain::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLawPlaneStrainTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestPowerLawPlaneStrain::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new PowerLawPlaneStrainTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestPowerLawPlaneStrain::testHasProperty(void)
{ // testHasProperty
  PowerLawPlaneStrain material;

  CPPUNIT_ASSERT(material.hasProperty("mu"));
  CPPUNIT_ASSERT(material.hasProperty("lambda"));
  CPPUNIT_ASSERT(material.hasProperty("density"));
  CPPUNIT_ASSERT(material.hasProperty("reference_strain_rate"));
  CPPUNIT_ASSERT(material.hasProperty("reference_stress"));
  CPPUNIT_ASSERT(material.hasProperty("power_law_exponent"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::materials::TestPowerLawPlaneStrain::testHasStateVar(void)
{ // testHasStateVar
  PowerLawPlaneStrain material;

  CPPUNIT_ASSERT(material.hasStateVar("stress_zz_initial"));
  CPPUNIT_ASSERT(material.hasStateVar("viscous_strain"));
  CPPUNIT_ASSERT(material.hasStateVar("stress4"));
  CPPUNIT_ASSERT(!material.hasStateVar("total_strain"));
} // testHasStateVar


// End of file 
