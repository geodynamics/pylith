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

#include "ObserverSolnStub.hh" // Implementation of class methods

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverSolnStub::ObserverSolnStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverSolnStub::~ObserverSolnStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify observer is compatible with solution.
void
pylith::problems::ObserverSolnStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    throw ObserverSolnStubException(ObserverSolnStubException::VERIFIED);
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Receive update (subject of observer).
void
pylith::problems::ObserverSolnStub::update(const PylithReal t,
                                           const PylithInt tindex,
                                           const pylith::topology::Field& solution,
                                           const bool infoOnly) {
    throw ObserverSolnStubException(ObserverSolnStubException::UPDATED);
} // update


// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverSolnStubException::ObserverSolnStubException(const MethodEnum value) :
    _methodCalled(value) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverSolnStubException::~ObserverSolnStubException(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Get method called.
pylith::problems::ObserverSolnStubException::MethodEnum
pylith::problems::ObserverSolnStubException::getMethodCalled(void) const {
    return _methodCalled;
} // getMethodCalled


// End of file
