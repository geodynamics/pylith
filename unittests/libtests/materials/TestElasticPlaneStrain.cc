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
  setupNormalizer();
} // setUp


// End of file 
