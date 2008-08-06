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
// Setup testing data.
void
pylith::materials::TestElasticStress1D::setUp(void)
{ // setUp
  _material = new ElasticStress1D();
  _matElastic = new ElasticStress1D();
  _data = new ElasticStress1DData();
  _dataElastic = new ElasticStress1DData();
} // setUp

// ----------------------------------------------------------------------
// Test usesUpdateProperties()
void
pylith::materials::TestElasticStress1D::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  ElasticStress1D material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test updateProperties()
void
pylith::materials::TestElasticStress1D::testUpdateProperties(void)
{ // testUpdateProperties
  ElasticStress1D material;
  ElasticStress1DData data;

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
    initialState[i] = 0.1*i;
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
