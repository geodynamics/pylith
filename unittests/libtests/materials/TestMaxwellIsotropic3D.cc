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

#include "TestMaxwellIsotropic3D.hh" // Implementation of class methods

#include "data/MaxwellIsotropic3DElasticData.hh" // USES MaxwellIsotropic3DElasticData
#include "data/MaxwellIsotropic3DTimeDepData.hh" // USES MaxwellIsotropic3DTimeDepData

#include "pylith/materials/MaxwellIsotropic3D.hh" // USES MaxwellIsotropic3D

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaxwellIsotropic3D::setUp(void)
{ // setUp
  _material = new MaxwellIsotropic3D();
  _matElastic = new MaxwellIsotropic3D();
  _data = new MaxwellIsotropic3DElasticData();
  _dataElastic = new MaxwellIsotropic3DElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  MaxwellIsotropic3D material;

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
pylith::materials::TestMaxwellIsotropic3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  MaxwellIsotropic3D material;

  material.useElasticBehavior(true);
  
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::MaxwellIsotropic3D::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic ==
		   material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestMaxwellIsotropic3D::testHasStateVars(void)
{ // testHasStateVars
  MaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
   _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestMaxwellIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestMaxwellIsotropic3D::testHasProperty(void)
{ // testHasProperty
  MaxwellIsotropic3D material;

  CPPUNIT_ASSERT(material.hasProperty("mu"));
  CPPUNIT_ASSERT(material.hasProperty("lambda"));
  CPPUNIT_ASSERT(material.hasProperty("density"));
  CPPUNIT_ASSERT(material.hasProperty("maxwell_time"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::materials::TestMaxwellIsotropic3D::testHasStateVar(void)
{ // testHasStateVar
  MaxwellIsotropic3D material;

  CPPUNIT_ASSERT(material.hasStateVar("total_strain"));
  CPPUNIT_ASSERT(material.hasStateVar("viscous_strain"));
  CPPUNIT_ASSERT(!material.hasStateVar("stress"));
} // testHasStateVar


// End of file 
