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

#include "TestElasticStrain1D.hh" // Implementation of class methods

#include "data/ElasticStrain1DData.hh" // USES ElasticStrain1DData

#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticStrain1D );

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestElasticStrain1D::testDBValues(void)
{ // testDBValues
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestElasticStrain1D::testParameters(void)
{ // testParameters
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestElasticStrain1D::testDBToParameters(void)
{ // testDBToParameters
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticStrain1D::testCalcDensity(void)
{ // testCalcDensity
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticStrain1D::testCalcStress(void)
{ // testCalcStress
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestElasticStrain1D::testCalcElasticConsts(void)
{ // testElasticConsts
  ElasticStrain1D material;
  ElasticStrain1DData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestElasticStrain1D::testUpdateState(void)
{ // testUpdateState
  ElasticStrain1D material;

  std::vector<double_array> totalStrain;
  material.updateState(totalStrain);
} // testUpdateState


// End of file 
