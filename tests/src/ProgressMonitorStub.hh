// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/problems/ProgressMonitor.hh" // ISA ProgressMonitor

class pylith::problems::ProgressMonitorStub : public pylith::problems::ProgressMonitor {
    friend class TestProgressMonitor; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ProgressMonitorStub(void);

    /// Destructor
    ~ProgressMonitorStub(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    void _open(void);

    /// Close progress monitor.
    void _close(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////

    struct State {
        double current; ///< Current time/step.
        time_t now; ///< Current time.
        double percentComplete; ///< Percent complete
        const char* finished; ///< Time stamp when finished
    };
    State _state; ///< Current state

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ProgressMonitorStub(const ProgressMonitorStub&); ///< Not implemented.
    const ProgressMonitorStub& operator=(const ProgressMonitorStub&); ///< Not implemented

}; // ProgressMonitorStub

// End of file
