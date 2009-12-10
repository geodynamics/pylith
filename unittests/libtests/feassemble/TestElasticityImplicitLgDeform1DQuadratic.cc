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

#include "TestElasticityImplicitLgDeform1DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityImplicitLgDeformData1DQuadratic.hh"

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D
#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform1DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform1DQuadratic::setUp(void)
{ // setUp
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData1DQuadratic();
  _gravityField = 0;
  GeometryLine1D geometry;
  CPPUNIT_ASSERT(0 != _quadrature);
  _quadrature->refGeometry(&geometry);
  _material = new materials::ElasticStrain1D;

  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticStrain1D"),
		       std::string(_data->matType));
} // setUp


// End of file 
