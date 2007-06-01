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

#include "TestElasticityImplicit2DLinear.hh" // Implementation of class methods

#include "data/ElasticityImplicitData2DLinear.hh"

#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature2D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit2DLinear );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit2DLinear::setUp(void)
{ // setUp
  _data = new ElasticityImplicitData2DLinear();
  _quadrature = new Quadrature2D();
  _material = new materials::ElasticPlaneStrain;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"),
		       std::string(_data->matType));
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityImplicit2DLinear::tearDown(void)
{ // tearDown
  delete _data;
  delete _quadrature;
  delete _material;
} // tearDown


// End of file 
