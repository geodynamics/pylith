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

#include "TestElasticityExplicitLgDeform1DQuadratic.hh" // Implementation of class methods

#include "data/ElasticityExplicitLgDeformData1DQuadratic.hh"

#include "pylith/topology/Mesh.hh" // USES Quadrature<Mesh>
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/GeometryLine1D.hh" // USES GeometryLine1D
#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestElasticityExplicitLgDeform1DQuadratic );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::feassemble::TestElasticityExplicitLgDeform1DQuadratic::setUp(void)
{ // setUp
  TestElasticityExplicitLgDeform::setUp();

  _data = new ElasticityExplicitLgDeformData1DQuadratic();
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
