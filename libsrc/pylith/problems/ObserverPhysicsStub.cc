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

#include "ObserverPhysicsStub.hh" // Implementation of class methods

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverPhysicsStub::ObserverPhysicsStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverPhysicsStub::~ObserverPhysicsStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify observer is compatible with solution.
void
pylith::problems::ObserverPhysicsStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    throw ObserverPhysicsStubException(ObserverPhysicsStubException::VERIFIED);
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Receive update (subject of observer).
void
pylith::problems::ObserverPhysicsStub::update(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution,
                                              const bool infoOnly) {
    throw ObserverPhysicsStubException(ObserverPhysicsStubException::UPDATED);
} // update


// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverPhysicsStubException::ObserverPhysicsStubException(const MethodEnum value) :
    _methodCalled(value) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverPhysicsStubException::~ObserverPhysicsStubException(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Get method called.
pylith::problems::ObserverPhysicsStubException::MethodEnum
pylith::problems::ObserverPhysicsStubException::getMethodCalled(void) const {
    return _methodCalled;
} // getMethodCalled


// End of file
