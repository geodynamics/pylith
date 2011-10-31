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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticityImplicitGrav2DLinear.hh" // Implementation of class methods

#include "data/ElasticityImplicitGravData2DLinear.hh"

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryTri2D.hh" // USES GeometryTri2D
#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityImplicitGrav2DLinear );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityImplicitGrav2DLinear::setUp(void)
{ // setUp
  TestElasticityImplicit::setUp();

  _data = new ElasticityImplicitGravData2DLinear();
  _gravityField = new spatialdata::spatialdb::GravityField();
  CPPUNIT_ASSERT(0 != _quadrature);
  CPPUNIT_ASSERT(0 != _gravityField);
  GeometryTri2D geometry;
  _quadrature->refGeometry(&geometry);

  const PylithScalar g = 1.0e8;
  const PylithScalar gravityDir[] = { 0.0, -1.0, 0.0 };
  _gravityField->gravAcceleration(g);
  _gravityField->gravityDir(gravityDir[0], gravityDir[1], gravityDir[2]);

  _material = new materials::ElasticPlaneStrain;
  CPPUNIT_ASSERT(0 != _material);
  
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT_EQUAL(std::string("ElasticPlaneStrain"),
		       std::string(_data->matType));
} // setUp


// End of file 
