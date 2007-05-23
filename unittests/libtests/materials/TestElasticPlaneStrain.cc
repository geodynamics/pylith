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

#include "TestElasticPlaneStrain.hh" // Implementation of class methods

#include "data/ElasticPlaneStrainData.hh" // USES ElasticPlaneStrainData

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticPlaneStrain );

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestElasticPlaneStrain::testDBValues(void)
{ // testDBValues
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestElasticPlaneStrain::testParameters(void)
{ // testParameters
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestElasticPlaneStrain::testDBToParameters(void)
{ // testDBToParameters
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticPlaneStrain::testCalcDensity(void)
{ // testCalcDensity
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticPlaneStrain::testCalcStress(void)
{ // testCalcStress
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestElasticPlaneStrain::testCalcElasticConsts(void)
{ // testElasticConsts
  ElasticPlaneStrain material;
  ElasticPlaneStrainData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestElasticPlaneStrain::testUpdateState(void)
{ // testUpdateState
  ElasticPlaneStrain material;

  std::vector<double_array> totalStrain;
  material.updateState(totalStrain);
} // testUpdateState


// End of file 
