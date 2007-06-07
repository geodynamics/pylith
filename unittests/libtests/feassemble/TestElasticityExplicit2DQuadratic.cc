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

#include "TestElasticityExplicit2DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityExplicitData2DQuadratic.hh"

#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature2D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit2DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit2DQuadratic::setUp(void)
{ // setUp
  _data = new ElasticityExplicitData2DQuadratic();
  _quadrature = new Quadrature2D();
  _material = new materials::ElasticPlaneStrain;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"),
		       std::string(_data->matType));
} // setUp


// End of file 
