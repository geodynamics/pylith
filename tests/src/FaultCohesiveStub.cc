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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "FaultCohesiveStub.hh" // implementation of object methods

#include "StubMethodTracker.hh" // USES StubMethodTracker
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveStub::FaultCohesiveStub(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveStub::~FaultCohesiveStub(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::verifyConfiguration");
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesiveStub::createIntegrator(const pylith::topology::Field& solution,
                                                    const std::vector<pylith::materials::Material*>& materials) {
    return NULL;
}

// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveStub::createAuxiliaryField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& physicsMesh) {
    pylith::testing::StubMethodTracker tracker("pylith::faults::FaultCohesiveStub::createAuxiliaryField");

    return NULL;
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesiveStub::_getAuxiliaryFactory(void) {
    return NULL;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::faults::FaultCohesiveStub::_setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                                       const pylith::topology::Field& solution,
                                                       const std::vector<pylith::materials::Material*>& materials) const {}


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::faults::FaultCohesiveStub::_setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                                       const pylith::topology::Field& solution,
                                                       const std::vector<pylith::materials::Material*>& materials) const {}



// End of file
