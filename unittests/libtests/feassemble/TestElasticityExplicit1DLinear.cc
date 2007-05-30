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

#include "TestElasticityExplicit1DLinear.hh" // Implementation of class methods

#include "data/ElasticityExplicitData1DLinear.hh"

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D
#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit1DLinear );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit1DLinear::setUp(void)
{ // setUp
  _data = new ElasticityExplicitData1DLinear();
  _quadrature = new Quadrature1D();
  _material = new materials::ElasticStrain1D;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"),
		       std::string(_data->matType));
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityExplicit1DLinear::tearDown(void)
{ // tearDown
  delete _data;
  delete _quadrature;
  delete _material;
} // tearDown


// End of file 
