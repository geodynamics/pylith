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

#include "TestElasticityImplicit3DLinear.hh" // Implementation of class methods

#include "data/ElasticityImplicitData3DLinear.hh"

#include "pylith/feassemble/Quadrature3D.hh" // USES Quadrature3D
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit3DLinear );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit3DLinear::setUp(void)
{ // setUp
  _data = new ElasticityImplicitData3DLinear();
  _quadrature = new Quadrature3D();
  _material = new materials::ElasticIsotropic3D;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"),
		       std::string(_data->matType));
} // setUp


// End of file 
