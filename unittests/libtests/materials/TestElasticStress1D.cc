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
  setupNormalizer();
} // setUp


// End of file 
