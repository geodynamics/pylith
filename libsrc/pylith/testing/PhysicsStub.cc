// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "PhysicsStub.hh" // Implementation of class methods

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::PhysicsStub::PhysicsStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::PhysicsStub::~PhysicsStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::PhysicsStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::verifyConfiguration");
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::problems::PhysicsStub::createIntegrator(const pylith::topology::Field& solution) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createIntegrator");

    return NULL;
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::problems::PhysicsStub::createConstraint(const pylith::topology::Field& solution) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createConstraint");

    return NULL;
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::problems::PhysicsStub::createAuxiliaryField(const pylith::topology::Field& solution,
                                                    const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createAuxiliaryField");

    return NULL;
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::problems::PhysicsStub::createDerivedField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::createDerivedField");

    return NULL;
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::problems::PhysicsStub::_getAuxiliaryFactory(void) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::PhysicsStub::_getAuxiliaryFactory");

    return NULL;
} // _getAuxiliaryFactory


// End of file
