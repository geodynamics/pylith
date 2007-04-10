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


// End of file 
