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

#include "ObserverPhysicsStub.hh" // Implementation of class methods

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverPhysicsStub::ObserverPhysicsStub(void) :
    _timeScale(1.0)
{}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ObserverPhysicsStub::~ObserverPhysicsStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set time scale.
void
pylith::problems::ObserverPhysicsStub::setTimeScale(const PylithReal value) {
    _timeScale = value;
} // setTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Get time scale.
PylithReal
pylith::problems::ObserverPhysicsStub::getTimeScale(void) const {
    return _timeScale;
} // getTimeScale


// ---------------------------------------------------------------------------------------------------------------------
// Verify observer is compatible with solution.
void
pylith::problems::ObserverPhysicsStub::verifyConfiguration(const pylith::topology::Field& solution) const {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ObserverPhysicsStub::verifyConfiguration");
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Receive update (subject of observer).
void
pylith::problems::ObserverPhysicsStub::update(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution,
                                              const bool infoOnly) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ObserverPhysicsStub::update");
} // update


// End of file
