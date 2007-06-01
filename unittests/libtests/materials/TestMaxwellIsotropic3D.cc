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

#include "data/MaxwellIsotropic3DData.hh" // USES MaxwellIsotropic3DData

#include "pylith/materials/MaxwellIsotropic3D.hh" // USES MaxwellIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestMaxwellIsotropic3D::testDBValues(void)
{ // testDBValues
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestMaxwellIsotropic3D::testParameters(void)
{ // testParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestMaxwellIsotropic3D::testDBToParameters(void)
{ // testDBToParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcDensity(void)
{ // testCalcDensity
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStress(void)
{ // testCalcStress
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConsts(void)
{ // testElasticConsts
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdateState(void)
{ // testUpdateState
  MaxwellIsotropic3D material;

  std::vector<double_array> totalStrain;
  material.updateState(totalStrain);
} // testUpdateState

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  MaxwellIsotropic3D material;

  CPPUNIT_ASSERT_EQUAL(false, material._needNewJacobian);

  const double dt1 = 1.0;
  material.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, material.timeStep());
  CPPUNIT_ASSERT_EQUAL(false, material.needNewJacobian());

  const double dt2 = 2.0;
  CPPUNIT_ASSERT_EQUAL(dt2, material.timeStep());
  CPPUNIT_ASSERT_EQUAL(true, material.needNewJacobian());
} // testTimeStep


// End of file 
