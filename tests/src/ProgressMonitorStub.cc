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

#include "ProgressMonitorStub.hh" // Implementation of class methods

#include "tests/src/StubMethodTracker.hh" // USES StubMethodTracker

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ProgressMonitorStub::ProgressMonitorStub(void) {
    _state.current = 0;
    _state.now = 0;
    _state.percentComplete = 0;
    _state.finished = NULL;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ProgressMonitorStub::~ProgressMonitorStub(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Open progress monitor.
void
pylith::problems::ProgressMonitorStub::_open(void) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ProgressMonitorStub::_open");
} // _open


// ------------------------------------------------------------------------------------------------
// Close progress monitor.
void
pylith::problems::ProgressMonitorStub::_close(void) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ProgressMonitorStub::_close");
} // _close


// End of file
