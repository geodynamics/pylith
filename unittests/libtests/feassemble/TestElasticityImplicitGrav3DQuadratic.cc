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

#include "TestElasticityImplicitGrav3DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityImplicitGravData3DQuadratic.hh"

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/Quadrature3D.hh" // USES Quadrature3D
#include "pylith/feassemble/GeometryTet3D.hh" // USES GeometryTet3D
#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav3DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav3DQuadratic::setUp(void)
{ // setUp
  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData3DQuadratic();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != _gravityField);
  GeometryTet3D geometry;
  _quadrature->refGeometry(&geometry);

  const double g = 1.0e8;
  _gravityField->gravAcceleration(g);

  _material = new materials::ElasticIsotropic3D;
  CPPUNIT_ASSERT(0 != _material);
  
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticIsotropic3D"),
		       std::string(_data->matType));
} // setUp


// End of file 
