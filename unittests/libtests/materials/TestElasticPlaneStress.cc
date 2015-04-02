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

#include "TestElasticPlaneStress.hh" // Implementation of class methods

#include "data/ElasticPlaneStressData.hh" // USES ElasticPlaneStressData

#include "pylith/materials/ElasticPlaneStress.hh" // USES ElasticPlaneStress

#include "pylith/utils/constdefs.h" // USES MAXSCALAR
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticPlaneStress );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticPlaneStress::setUp(void)
{ // setUp
  _material = new ElasticPlaneStress();
  _matElastic = new ElasticPlaneStress();
  _data = new ElasticPlaneStressData();
  _dataElastic = new ElasticPlaneStressData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test stableTimeStepImplicit().
void
pylith::materials::TestElasticPlaneStress::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  assert(0 != _matElastic);

  topology::Mesh mesh;

  const PylithScalar dt = _matElastic->stableTimeStepImplicit(mesh);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), dt/pylith::PYLITH_MAXSCALAR);
} // testStableTimeStepImplicit


// End of file 
