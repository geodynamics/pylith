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
// Setup testing data.
void
pylith::materials::TestElasticPlaneStrain::setUp(void)
{ // setUp
  _material = new ElasticPlaneStrain();
  _matElastic = new ElasticPlaneStrain();
  _data = new ElasticPlaneStrainData();
  _dataElastic = new ElasticPlaneStrainData();
} // setUp

// ----------------------------------------------------------------------
// Test usesUpdateProperties()
void
pylith::materials::TestElasticPlaneStrain::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  ElasticPlaneStrain material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test updateProperties()
void
pylith::materials::TestElasticPlaneStrain::testUpdateProperties(void)
{ // testUpdateProperties
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
  const int initialStateSize = 9;
  double_array totalStrain(tensorSize);
  double_array initialState(initialStateSize);
  for (int i=0; i < tensorSize; ++i) {
    totalStrain[i] = i;
    initialState[i] = 0;
  } // for
  
  material._updateProperties(&parameters[0], numParams, &totalStrain[0],
		  tensorSize, &initialState[0], initialStateSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParams; ++i)
    for (int j=0; j < paramsSize; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i+j), parameters[i*paramsSize+j],
				   tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
} // testUpdateProperties


// End of file 
