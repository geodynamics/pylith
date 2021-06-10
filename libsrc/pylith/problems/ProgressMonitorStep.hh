// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2019 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/*** @file libsrc/problems/ProgressMonitorStep.hh
 *
 * @brief Abstract base class for objects defining physics, such as behavior
 * of a bulk makterial, boundary condition, interface, or constraint.
 */

#if !defined(pylith_problems_ProgressMonitorStep_hh)
#define pylith_problems_ProgressMonitorStep_hh

#include "ProgressMonitor.hh" // ISA ProgressMonitor

#include <fstream> // HASA std::ofstream

class pylith::problems::ProgressMonitorStep : public pylith::problems::ProgressMonitor {
    friend class TestProgressMonitorStep; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    ProgressMonitorStep(void);

    /// Destructor
    virtual ~ProgressMonitorStep(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    void _open(void);

    /// Close progress monitor.
    void _close(void);

    /** Update progress.
     *
     * @param[in] current step.
     * @param[in] now Current time.
     * @param[in] percentComplete Percent completed
     * @param[in] finished Time stamp of estimated finish.
     */
    void _update(const int current,
                 const time_t& now,
                 const double percentComplete,
                 const char* finished);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    std::ofstream _sout; ///< Output stream.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ProgressMonitorStep(const ProgressMonitorStep&); ///< Not implemented.
    const ProgressMonitorStep& operator=(const ProgressMonitorStep&); ///< Not implemented.

}; // ProgressMonitorStep

#endif // pylith_problems_ProgressMonitorStep_hh

// End of file
