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

#include "TestElasticityExplicit3DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityExplicitData3DQuadratic.hh"

#include "pylith/feassemble/Quadrature3D.hh" // USES Quadrature3D
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit3DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit3DQuadratic::setUp(void)
{ // setUp
  _data = new ElasticityExplicitData3DQuadratic();
  _quadrature = new Quadrature3D();
  _material = new materials::ElasticIsotropic3D;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"),
		       std::string(_data->matType));
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::feassemble::TestElasticityExplicit3DQuadratic::tearDown(void)
{ // tearDown
  delete _data;
  delete _quadrature;
  delete _material;
} // tearDown


// End of file 
