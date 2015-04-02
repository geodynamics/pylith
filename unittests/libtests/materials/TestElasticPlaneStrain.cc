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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticPlaneStrain.hh" // Implementation of class methods

#include "data/ElasticPlaneStrainData.hh" // USES ElasticPlaneStrainData

#include "pylith/materials/ElasticPlaneStrain.hh" // USES ElasticPlaneStrain

#include "pylith/utils/constdefs.h" // USES MAXSCALAR
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticPlaneStrain );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticPlaneStrain::setUp(void)
{ // setUp
  _material = new ElasticPlaneStrain();
  _matElastic = new ElasticPlaneStrain();
  _data = new ElasticPlaneStrainData();
  _dataElastic = new ElasticPlaneStrainData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test stableTimeStepImplicit().
void
pylith::materials::TestElasticPlaneStrain::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  assert(0 != _matElastic);

  topology::Mesh mesh;

  const PylithScalar dt = _matElastic->stableTimeStepImplicit(mesh);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), dt/pylith::PYLITH_MAXSCALAR);
} // testStableTimeStepImplicit


// End of file 
