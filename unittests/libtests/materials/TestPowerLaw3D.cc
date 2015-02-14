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

#include "TestPowerLaw3D.hh" // Implementation of class methods

#include "data/PowerLaw3DElasticData.hh" // USES PowerLaw3DElasticData
#include "data/PowerLaw3DTimeDepData.hh" // USES PowerLaw3DTimeDepData

#include "pylith/materials/PowerLaw3D.hh" // USES PowerLaw3D

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestPowerLaw3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestPowerLaw3D::setUp(void)
{ // setUp
  _material = new PowerLaw3D();
  _matElastic = new PowerLaw3D();
  _data = new PowerLaw3DElasticData();
  _dataElastic = new PowerLaw3DElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestPowerLaw3D::testTimeStep(void)
{ // testTimeStep
  PowerLaw3D material;

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
pylith::materials::TestPowerLaw3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  PowerLaw3D material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::PowerLaw3D::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::PowerLaw3D::_updateStateVarsViscoelastic ==
		   material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestPowerLaw3D::testHasStateVars(void)
{ // testHasStateVars
  PowerLaw3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestPowerLaw3D::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestPowerLaw3D::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestPowerLaw3D::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestPowerLaw3D::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLaw3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestPowerLaw3D::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLaw3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestPowerLaw3D::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new PowerLaw3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestPowerLaw3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new PowerLaw3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestPowerLaw3D::testHasProperty(void)
{ // testHasProperty
  PowerLaw3D material;

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
pylith::materials::TestPowerLaw3D::testHasStateVar(void)
{ // testHasStateVar
  PowerLaw3D material;

  CPPUNIT_ASSERT(material.hasStateVar("viscous_strain"));
  CPPUNIT_ASSERT(material.hasStateVar("stress"));
  CPPUNIT_ASSERT(!material.hasStateVar("total_strain"));
} // testHasStateVar


// End of file 
