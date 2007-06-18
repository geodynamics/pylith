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

#include "TestElasticIsotropic3D.hh" // Implementation of class methods

#include "data/ElasticIsotropic3DData.hh" // USES ElasticIsotropic3DData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticIsotropic3D );

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestElasticIsotropic3D::testUsesUpdateState(void)
{ // testUsesUpdateState
  ElasticIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateState());
} // testUsesUpdateState

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestElasticIsotropic3D::testDBValues(void)
{ // testDBValues
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestElasticIsotropic3D::testParameters(void)
{ // testParameters
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestElasticIsotropic3D::testDBToParameters(void)
{ // testDBToParameters
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestElasticIsotropic3D::testCalcDensity(void)
{ // testCalcDensity
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStress()
void
pylith::materials::TestElasticIsotropic3D::testCalcStress(void)
{ // testCalcStress
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testCalcStress(&material, data);
} // testCalcStress

// ----------------------------------------------------------------------
// Test calcElasticConsts()
void
pylith::materials::TestElasticIsotropic3D::testCalcElasticConsts(void)
{ // testElasticConsts
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;
  _testCalcElasticConsts(&material, data);
} // testElasticConsts

// ----------------------------------------------------------------------
// Test updateState()
void
pylith::materials::TestElasticIsotropic3D::testUpdateState(void)
{ // testUpdateState
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;

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
