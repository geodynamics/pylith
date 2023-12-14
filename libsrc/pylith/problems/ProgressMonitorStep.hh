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

#include "pylith/problems/ProgressMonitor.hh" // ISA ProgressMonitor

#include <fstream> // HASA std::ofstream

class pylith::problems::ProgressMonitorStep : public pylith::problems::ProgressMonitor {
    friend class TestProgressMonitorStep; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    ProgressMonitorStep(void);

    /// Destructor
    virtual ~ProgressMonitorStep(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Update progress.
     *
     * @param[in] current Current step.
     * @param[in] start Starting step.
     * @param[in] stop Ending step.
     */
    void update(const size_t current,
                const size_t start,
                const size_t stop);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    void _open(void);

    /// Close progress monitor.
    void _close(void);

    /** Update progress.
     *
     * @param[in] setp Current step.
     * @param[in] now Current date/time.
     * @param[in] percentComplete Percent completed
     * @param[in] finished Time stamp of estimated finish.
     */
    void _update(const size_t step,
                 const time_t& now,
                 const double percentComplete,
                 const char* finished);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::ofstream _sout; ///< Output stream.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ProgressMonitorStep(const ProgressMonitorStep&); ///< Not implemented.
    const ProgressMonitorStep& operator=(const ProgressMonitorStep&); ///< Not implemented.

}; // ProgressMonitorStep

// End of file
