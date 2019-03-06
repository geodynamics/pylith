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

#include "FaultCohesiveStub.hh" // implementation of object methods

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveStub::FaultCohesiveStub(void)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveStub::~FaultCohesiveStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    throw FaultCohesiveStubException(FaultCohesiveStubException::VERIFY_CONFIGURATION);
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesiveStub::createIntegrator(const pylith::topology::Field& solution) {
    throw FaultCohesiveStubException(FaultCohesiveStubException::CREATE_INTEGRATOR);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::faults::FaultCohesiveStub::createConstraint(const pylith::topology::Field& solution) {
    throw FaultCohesiveStubException(FaultCohesiveStubException::CREATE_CONSTRAINT);
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesiveStub::createAuxiliaryField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& physicsMesh) {
    throw FaultCohesiveStubException(FaultCohesiveStubException::CREATE_AUXILIARY_FIELD);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::faults::FaultCohesiveStub::createDerivedField(const pylith::topology::Field& solution,
                                                      const pylith::topology::Mesh& physicsMesh) {
    throw FaultCohesiveStubException(FaultCohesiveStubException::CREATE_DERIVED_FIELD);
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


// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::faults::FaultCohesiveStubException::FaultCohesiveStubException(const MethodEnum value) :
    _methodCalled(value) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::faults::FaultCohesiveStubException::~FaultCohesiveStubException(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Get method called.
pylith::faults::FaultCohesiveStubException::MethodEnum
pylith::faults::FaultCohesiveStubException::getMethodCalled(void) const {
    return _methodCalled;
} // getMethodCalled


// End of file
