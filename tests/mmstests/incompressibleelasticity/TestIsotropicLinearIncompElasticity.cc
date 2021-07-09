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

#include "TestIsotropicLinearIncompElasticity.hh" // Implementation of class methods

#include "pylith/materials/IncompressibleElasticity.hh" // USES IncompressibleElasticity
#include "pylith/materials/IsotropicLinearIncompElasticity.hh" // USES IsotropicLinearIncompElasticity

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestIsotropicLinearIncompElasticity::setUp(void) {
    TestIncompressibleElasticity::setUp();

    _rheology = new pylith::materials::IsotropicLinearIncompElasticity();CPPUNIT_ASSERT(_rheology);

    CPPUNIT_ASSERT(_material);
    _material->setBulkRheology(_rheology);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestIsotropicLinearIncompElasticity::tearDown(void) {
    delete _rheology;_rheology = NULL;

    TestIncompressibleElasticity::tearDown();
} // tearDown


// End of file
