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
// Test usesUpdateState()
void
pylith::materials::TestElasticPlaneStrain::testUsesUpdateState(void)
{ // testUsesUpdateState
  ElasticPlaneStrain material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateState());
} // testUsesUpdateState

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
  ElasticPlaneStrainData data;

  const int numParams = data.numParameters;

  double_array parameters(numParams);
  const int paramsSize = 1;
  for (int i=0; i < numParams; ++i) {
    for (int j=0; j < paramsSize; ++j)
      parameters[i*paramsSize+j] = i+j;
  } // for
    
  const int tensorSize = 9;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;
  
  material._updateState(&parameters[0], numParams, &totalStrain[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParams; ++i)
    for (int j=0; j < paramsSize; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i+j), parameters[i*paramsSize+j],
				   tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
} // testUpdateState


// End of file 
