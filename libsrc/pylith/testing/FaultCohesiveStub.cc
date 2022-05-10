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

#include "FaultCohesiveStub.hh" // implementation of object methods

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveStub::FaultCohesiveStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveStub::~FaultCohesiveStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::verifyConfiguration");
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesiveStub::createIntegrator(const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::createIntegrator");

    return NULL;
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::faults::FaultCohesiveStub::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::createConstraints");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveStub::createAuxiliaryField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::createAuxiliaryField");

    return NULL;
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::faults::FaultCohesiveStub::createDerivedField(const pylith::topology::Field& solution,
                                                      const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::createDerivedField");

    return NULL;
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveStub::_getAuxiliaryFactory(void) {
    return NULL;
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::FaultCohesiveStub::_updateKernelConstants(const PylithReal dt) {}


// End of file
