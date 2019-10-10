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

#include "ProgressMonitorStub.hh" // Implementation of class methods

#include "pylith/testing/StubMethodTracker.hh" // USES StubMethodTracker

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::problems::ProgressMonitorStub::ProgressMonitorStub(void) {
    _state.current = 0;
    _state.now = NULL;
    _state.percentComplete = 0;
    _state.finished = NULL;
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::ProgressMonitorStub::~ProgressMonitorStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Open progress monitor.
void
pylith::problems::ProgressMonitorStub::_open(void) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ProgressMonitorStub::_open");
} // _open


// ---------------------------------------------------------------------------------------------------------------------
// Close progress monitor.
void
pylith::problems::ProgressMonitorStub::_close(void) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ProgressMonitorStub::_close");
} // _close


// ---------------------------------------------------------------------------------------------------------------------
// Update progress.
void
pylith::problems::ProgressMonitorStub::_update(const double current,
                                               const time_t& now,
                                               const double percentComplete,
                                               const char* finished) {
    pylith::testing::StubMethodTracker tracker("pylith::problems::ProgressMonitorStub::_update");
    _state.current = current;
    _state.now = now;
    _state.percentComplete = percentComplete;
    _state.finished = finished;
} // _update


// End of file
