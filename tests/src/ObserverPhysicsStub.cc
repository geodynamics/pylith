// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "ObserverPhysicsStub.hh" // Implementation of class methods

#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ObserverPhysicsStub::ObserverPhysicsStub(void) :
    _timeScale(1.0) {}


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
                                              const NotificationType notification) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ObserverPhysicsStub::update");
} // update


// End of file
