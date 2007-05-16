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

#include "TestElasticPlaneStress.hh" // Implementation of class methods

#include "data/ElasticPlaneStressData.hh" // USES ElasticPlaneStressData

#include "pylith/materials/ElasticPlaneStress.hh" // USES ElasticPlaneStress

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticPlaneStress );

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestElasticPlaneStress::testDBValues(void)
{ // testDBValues
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestElasticPlaneStress::testParameters(void)
{ // testParameters
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestElasticPlaneStress::testDBToParameters(void)
{ // testDBToParameters
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticPlaneStress::testCalcDensity(void)
{ // testCalcDensity
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticPlaneStress::testCalcStress(void)
{ // testCalcStress
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestElasticPlaneStress::testCalcElasticConsts(void)
{ // testElasticConsts
  ElasticPlaneStress material;
  ElasticPlaneStressData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestElasticPlaneStress::testUpdateState(void)
{ // testUpdateState
  ElasticPlaneStress material;
  material.updateState();
} // testUpdateState


// End of file 
