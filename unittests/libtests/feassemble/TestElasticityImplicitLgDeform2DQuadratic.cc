// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticityImplicitLgDeform2DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityImplicitLgDeformData2DQuadratic.hh"

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitLgDeform2DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitLgDeform2DQuadratic::setUp(void)
{ // setUp
  TestElasticityImplicitLgDeform::setUp();

  _data = new ElasticityImplicitLgDeformData2DQuadratic();
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
