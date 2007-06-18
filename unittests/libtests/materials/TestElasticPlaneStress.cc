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
// Test usesUpdateState()
void
pylith::materials::TestElasticPlaneStress::testUsesUpdateState(void)
{ // testUsesUpdateState
  ElasticPlaneStress material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateState());
} // testUsesUpdateState

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
  ElasticPlaneStressData data;

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
