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
  setupNormalizer();
} // setUp


// End of file 
