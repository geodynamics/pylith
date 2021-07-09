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

/**
 * @file libsrc/problems/ProgressMonitorStub.hh
 *
 * @brief Minimal C++ implementation of ProgressMonitor to allow testing of basic ProgressMonitor functionality and use
 * of ProgressMoitor objects in other tests.
 */

#if !defined(pylith_problems_progressmonitorstub_hh)
#define pylith_problems_progressmonitorstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/problems/ProgressMonitor.hh" // ISA ProgressMonitor

class pylith::problems::ProgressMonitorStub : public pylith::problems::ProgressMonitor {
    friend class TestProgressMonitor; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ProgressMonitorStub(void);

    /// Destructor
    ~ProgressMonitorStub(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    void _open(void);

    /// Close progress monitor.
    void _close(void);

    /** Update progress.
     *
     * @param[in current Current step.
     * @param[in] now Current time.
     * @param[in] percentComplete Percent completed
     * @param[in] finished Time stamp of estimated finish.
     */
    void _update(const double current,
                 const time_t& now,
                 const double percentComplete,
                 const char* finished);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////

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

#endif // pylith_problems_progressmonitorstub_hh

// End of file
