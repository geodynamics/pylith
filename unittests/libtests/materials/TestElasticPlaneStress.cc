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
  setupNormalizer();
} // setUp


// End of file 
