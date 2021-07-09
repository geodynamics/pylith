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

/*** @file libsrc/problems/ProgressMonitorTime.hh
 *
 * @brief Abstract base class for objects defining physics, such as behavior
 * of a bulk makterial, boundary condition, interface, or constraint.
 */

#if !defined(pylith_problems_progressmonitortime_hh)
#define pylith_problems_progressmonitortime_hh

#include "ProgressMonitor.hh" // ISA ProgressMonitor

#include <fstream> // HASA std::ofstream

class pylith::problems::ProgressMonitorTime : public pylith::problems::ProgressMonitor {
    friend class TestProgressMonitorTime; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    ProgressMonitorTime(void);

    /// Destructor
    virtual ~ProgressMonitorTime(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set unit for simulation time in output.
     *
     * @param[in] Unit of time.
     */
    void setTimeUnit(const char* value);

    /** Set unit for simulation time in output.
     *
     * @param[in] Unit of time.
     */
    const char* getTimeUnit(void) const;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
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

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    double _baseTime; ///< Units of time as seconds.
    std::string _baseUnit; ///< Unit of time.
    std::ofstream _sout; ///< Output stream.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ProgressMonitorTime(const ProgressMonitorTime&); ///< Not implemented.
    const ProgressMonitorTime& operator=(const ProgressMonitorTime&); ///< Not implemented.

}; // ProgressMonitorTime

#endif // pylith_problems_progressmonitortime_hh

// End of file
