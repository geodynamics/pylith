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

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicit2DLinear );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicit2DLinear::setUp(void)
{ // setUp
  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitData2DLinear();
  _gravityField = 0;
  CPPUNIT_ASSERT(0 != _quadrature);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(0 != _material);

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"),
		       std::string(_data->matType));
} // setUp


// End of file 
