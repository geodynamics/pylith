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

#include "ObserverStub.hh" // Implementation of class methods

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::feassemble::ObserverStub::ObserverStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::ObserverStub::~ObserverStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Verify observer is compatible with solution.
void
pylith::feassemble::ObserverStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    throw ObserverStubException(ObserverStubException::VERIFIED);
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Receive update (subject of observer).
void
pylith::feassemble::ObserverStub::update(const PylithReal t,
                                         const PylithInt tindex,
                                         const pylith::topology::Field& solution,
                                         const bool infoOnly) {
    throw ObserverStubException(ObserverStubException::UPDATED);
} // update


// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::feassemble::ObserverStubException::ObserverStubException(const MethodEnum value) :
    _methodCalled(value) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::feassemble::ObserverStubException::~ObserverStubException(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Get method called.
pylith::feassemble::ObserverStubException::MethodEnum
pylith::feassemble::ObserverStubException::getMethodCalled(void) const {
    return _methodCalled;
} // getMethodCalled


// End of file
