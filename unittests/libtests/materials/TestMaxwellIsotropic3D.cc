// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMaxwellIsotropic3D.hh" // Implementation of class methods

#include "data/MaxwellIsotropic3DElasticData.hh" // USES MaxwellIsotropic3DElasticData
#include "data/MaxwellIsotropic3DTimeDepData.hh" // USES MaxwellIsotropic3DTimeDepData

#include "pylith/materials/MaxwellIsotropic3D.hh" // USES MaxwellIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  MaxwellIsotropic3D material;

  CPPUNIT_ASSERT_EQUAL(false, material._needNewJacobian);

  const double dt1 = 1.0;
  material.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(false, material.needNewJacobian());

  const double dt2 = 2.0;
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
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressElastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateElastic,
		       material._updateStateFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateViscoelastic,
		       material._updateStateFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestMaxwellIsotropic3D::testUsesUpdateState(void)
{ // testUsesUpdateState
  MaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.usesUpdateState());
} // testUsesUpdateState

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestMaxwellIsotropic3D::testDBValues(void)
{ // testDBValues
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestMaxwellIsotropic3D::testParameters(void)
{ // testParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestMaxwellIsotropic3D::testDBToParameters(void)
{ // testDBToParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcDensity(void)
{ // testCalcDensity
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressElastic(void)
{ // testCalcStressElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testCalcStress(&material, data);
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsElastic(void)
{ // testElasticConstsElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdateStateElastic(void)
{ // testUpdateStateElastic
  MaxwellIsotropic3D material;

  std::vector<double_array> totalStrain;
  bool elasFlag = true;
  material.useElasticBehavior(elasFlag);
  material.updateState(totalStrain);
} // testUpdateStateElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressTimeDep(void)
{ // testCalcStressTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;
  bool elasFlag = false;
  material.useElasticBehavior(elasFlag);
  double dt = 2.0e5;
  material.timeStep(dt);
  _testCalcStress(&material, data);
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;
  bool elasFlag = false;
  material.useElasticBehavior(elasFlag);
  double dt = 2.0e5;
  material.timeStep(dt);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdateStateTimeDep(void)
{ // testUpdateStateTimeDep
  MaxwellIsotropic3D material;

  std::vector<double_array> totalStrain;
  bool elasFlag = false;
  material.useElasticBehavior(elasFlag);
  double dt = 2.0e5;
  material.timeStep(dt);
  material.updateState(totalStrain);
} // testUpdateStateTimeDep


// End of file 
