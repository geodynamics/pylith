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

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicit3DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicit3DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicit::setUp();

  _data = new ElasticityExplicitData3DQuadratic();
  CPPUNIT_ASSERT(0 != _quadrature);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(0 != _material);
  
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"),
		       std::string(_data->matType));
} // setUp


// End of file 
