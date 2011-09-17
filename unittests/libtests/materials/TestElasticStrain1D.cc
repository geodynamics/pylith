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

#include "TestElasticStrain1D.hh" // Implementation of class methods

#include "data/ElasticStrain1DData.hh" // USES ElasticStrain1DData

#include "pylith/materials/ElasticStrain1D.hh" // USES ElasticStrain1D

#include "pylith/utils/constdefs.h" // USES MAXSCALAR
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticStrain1D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticStrain1D::setUp(void)
{ // setUp
  _material = new ElasticStrain1D();
  _matElastic = new ElasticStrain1D();
  _data = new ElasticStrain1DData();
  _dataElastic = new ElasticStrain1DData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test stableTimeStepImplicit().
void
pylith::materials::TestElasticStrain1D::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  assert(0 != _matElastic);

  topology::Mesh mesh;

  const PylithScalar dt = _matElastic->stableTimeStepImplicit(mesh);
  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0), dt/pylith::PYLITH_MAXSCALAR);
} // testStableTimeStepImplicit


// End of file 
