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
// Setup testing data.
void
pylith::materials::TestElasticPlaneStress::setUp(void)
{ // setUp
  _material = new ElasticPlaneStress();
  _matElastic = new ElasticPlaneStress();
  _data = new ElasticPlaneStressData();
  _dataElastic = new ElasticPlaneStressData();
} // setUp

// ----------------------------------------------------------------------
// Test usesUpdateProperties()
void
pylith::materials::TestElasticPlaneStress::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  ElasticPlaneStress material;
  CPPUNIT_ASSERT_EQUAL(false, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test updateProperties()
void
pylith::materials::TestElasticPlaneStress::testUpdateProperties(void)
{ // testUpdateProperties
  ElasticPlaneStress material;
  ElasticPlaneStressData data;

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
  
  material._updateProperties(&parameters[0], numParams, &totalStrain[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParams; ++i)
    for (int j=0; j < paramsSize; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i+j), parameters[i*paramsSize+j],
				   tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
} // testUpdateProperties


// End of file 
