// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearElasticity.hh" // Implementation of class methods

#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestIsotropicLinearElasticity::setUp(void) {
    TestElasticity::setUp();

    _rheology = new pylith::materials::IsotropicLinearElasticity();CPPUNIT_ASSERT(_rheology);

    CPPUNIT_ASSERT(_material);
    _material->setBulkRheology(_rheology);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestIsotropicLinearElasticity::tearDown(void) {
    delete _rheology;_rheology = NULL;

    TestElasticity::tearDown();
} // tearDown


// End of file
