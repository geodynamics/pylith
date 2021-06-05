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

#include "TestIsotropicLinearPoroelasticity.hh" // Implementation of class methods

#include "pylith/materials/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestIsotropicLinearPoroelasticity::setUp(void) {
    TestPoroelasticity::setUp();

    _rheology = new pylith::materials::IsotropicLinearPoroelasticity();CPPUNIT_ASSERT(_rheology);

    CPPUNIT_ASSERT(_material);
    _material->setBulkRheology(_rheology);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestIsotropicLinearPoroelasticity::tearDown(void) {
    delete _rheology;_rheology = NULL;

    TestPoroelasticity::tearDown();
} // tearDown


// End of file
