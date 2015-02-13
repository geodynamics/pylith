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

#include "TestDruckerPrager3D.hh" // Implementation of class methods

#include "data/DruckerPrager3DElasticData.hh" // USES DruckerPrager3DElasticData
#include "data/DruckerPrager3DTimeDepData.hh" // USES DruckerPrager3DTimeDepData

#include "pylith/materials/DruckerPrager3D.hh" // USES DruckerPrager3D

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestDruckerPrager3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestDruckerPrager3D::setUp(void)
{ // setUp
  _material = new DruckerPrager3D();
  DruckerPrager3D* matTmp = new DruckerPrager3D();
  CPPUNIT_ASSERT(matTmp);
  matTmp->allowTensileYield(true);
  _matElastic = matTmp;

  _data = new DruckerPrager3DElasticData();
  _dataElastic = new DruckerPrager3DElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestDruckerPrager3D::testTimeStep(void)
{ // testTimeStep
  DruckerPrager3D material;

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
pylith::materials::TestDruckerPrager3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  DruckerPrager3D material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::DruckerPrager3D::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_calcStressElastoplastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_calcElasticConstsElastoplastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::DruckerPrager3D::_updateStateVarsElastoplastic ==
		   material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test allowTensileYield()
void
pylith::materials::TestDruckerPrager3D::testAllowTensileYield(void)
{ // testAllowTensileYield
  DruckerPrager3D material;
  CPPUNIT_ASSERT(!material._allowTensileYield);

  material.allowTensileYield(true);
  CPPUNIT_ASSERT(material._allowTensileYield);
} // testAllowTensileYield

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestDruckerPrager3D::testHasStateVars(void)
{ // testHasStateVars
  DruckerPrager3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestDruckerPrager3D::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestDruckerPrager3D::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestDruckerPrager3D::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_updateStateVars();
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestDruckerPrager3D::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPrager3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestDruckerPrager3D::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPrager3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestDruckerPrager3D::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new DruckerPrager3DTimeDepData();

  PylithScalar dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_updateStateVars();

} // test_updateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestDruckerPrager3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new DruckerPrager3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestDruckerPrager3D::testHasProperty(void)
{ // testHasProperty
  DruckerPrager3D material;

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
pylith::materials::TestDruckerPrager3D::testHasStateVar(void)
{ // testHasStateVar
  DruckerPrager3D material;

  CPPUNIT_ASSERT(material.hasStateVar("plastic_strain"));
} // testHasStateVar

// ----------------------------------------------------------------------
// Test _dbToProperties().
void
pylith::materials::TestDruckerPrager3D::testDBToProperties(void)
{ // testDBToProperties
  CPPUNIT_ASSERT(0 != _material);
  
  TestMaterial::testDBToProperties();

  CPPUNIT_ASSERT_EQUAL(false, _material->isJacobianSymmetric());
} // testDBToProperties


// End of file 
