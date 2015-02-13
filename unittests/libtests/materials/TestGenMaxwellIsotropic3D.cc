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

#include "TestGenMaxwellIsotropic3D.hh" // Implementation of class methods

#include "data/GenMaxwellIsotropic3DElasticData.hh" // USES GenMaxwellIsotropic3DElasticData
#include "data/GenMaxwellIsotropic3DTimeDepData.hh" // USES GenMaxwellIsotropic3DTimeDepData

#include "pylith/materials/GenMaxwellIsotropic3D.hh" // USES GenMaxwellIsotropic3D

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestGenMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestGenMaxwellIsotropic3D::setUp(void)
{ // setUp
  _material = new GenMaxwellIsotropic3D();
  _matElastic = new GenMaxwellIsotropic3D();

  _data = new GenMaxwellIsotropic3DElasticData();
  _dataElastic = new GenMaxwellIsotropic3DElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  GenMaxwellIsotropic3D material;

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
pylith::materials::TestGenMaxwellIsotropic3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  GenMaxwellIsotropic3D material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic == 
		   material._calcElasticConstsFn);
  CPPUNIT_ASSERT(&pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsViscoelastic ==
		 material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestGenMaxwellIsotropic3D::testHasStateVars(void)
{ // testHasStateVars
  GenMaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcStressElastic(void)
{ // testCalcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcElasticConstsElastic(void)
{ // testElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateVarsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_updateStateVarsElastic(void)
{ // testUpdateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // testUpdateStateVarsElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcStressTimeDep(void)
{ // testCalcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateVarsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_updateStateVarsTimeDep(void)
{ // testUpdateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // testUpdateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit


// End of file 
