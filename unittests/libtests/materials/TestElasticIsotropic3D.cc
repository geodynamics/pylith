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
// Setup testing data.
void
pylith::materials::TestElasticIsotropic3D::setUp(void)
{ // setUp
  _material = new ElasticIsotropic3D();
  _matElastic = new ElasticIsotropic3D();
  _data = new ElasticIsotropic3DData();
  _dataElastic = new ElasticIsotropic3DData();
} // setUp

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestElasticIsotropic3D::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  ElasticIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test updateProperties()
void
pylith::materials::TestElasticIsotropic3D::testUpdateProperties(void)
{ // testUpdateProperties
  ElasticIsotropic3D material;
  ElasticIsotropic3DData data;

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
