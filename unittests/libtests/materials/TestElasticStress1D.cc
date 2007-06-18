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

#include "TestElasticStress1D.hh" // Implementation of class methods

#include "data/ElasticStress1DData.hh" // USES ElasticStress1DData

#include "pylith/materials/ElasticStress1D.hh" // USES ElasticStress1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticStress1D );

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestElasticStress1D::testUsesUpdateState(void)
{ // testUsesUpdateState
  ElasticStress1D material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateState());
} // testUsesUpdateState

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestElasticStress1D::testDBValues(void)
{ // testDBValues
  ElasticStress1D material;
  ElasticStress1DData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestElasticStress1D::testParameters(void)
{ // testParameters
  ElasticStress1D material;
  ElasticStress1DData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestElasticStress1D::testDBToParameters(void)
{ // testDBToParameters
  ElasticStress1D material;
  ElasticStress1DData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticStress1D::testCalcDensity(void)
{ // testCalcDensity
  ElasticStress1D material;
  ElasticStress1DData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticStress1D::testCalcStress(void)
{ // testCalcStress
  ElasticStress1D material;
  ElasticStress1DData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestElasticStress1D::testCalcElasticConsts(void)
{ // testElasticConsts
  ElasticStress1D material;
  ElasticStress1DData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestElasticStress1D::testUpdateState(void)
{ // testUpdateState
  ElasticStress1D material;
  ElasticStress1DData data;

  const int numParams = data.numParameters;

  std::vector<double_array> parameters(numParams);
  const int paramsSize = 1;
  for (int i=0; i < numParams; ++i) {
    parameters[i].resize(numParams);
    for (int j=0; j < paramsSize; ++j)
      parameters[i][j] = i+j;
  } // for
    
  const int tensorSize = 9;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;
  
  material._updateState(&parameters, totalStrain);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParams; ++i)
    for (int j=0; j < paramsSize; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i+j), parameters[i][j], tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
} // testUpdateState


// End of file 
