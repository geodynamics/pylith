// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "PhysicsStub.hh" // Implementation of class methods

#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::PhysicsStub::PhysicsStub(void) :
    _auxiliaryFactory(new pylith::feassemble::AuxiliaryFactory) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::PhysicsStub::~PhysicsStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::PhysicsStub::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    Physics::deallocate();

    PYLITH_METHOD_END;
} // deallocate


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
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// End of file
