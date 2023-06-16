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
