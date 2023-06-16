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

/*** @file libsrc/problems/ProgressMonitor.hh
 *
 * @brief Abstract base class for objects defining physics, such as behavior
 * of a bulk makterial, boundary condition, interface, or constraint.
 */

#if !defined(pylith_problems_progressmonitor_hh)
#define pylith_problems_progressmonitor_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include <string> // HASA std::string
#include <ctime> // HASA time_t

class pylith::problems::ProgressMonitor : public pylith::utils::PyreComponent {
    friend class TestProgressMonitor; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    ProgressMonitor(void);

    /// Destructor
    virtual ~ProgressMonitor(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set how often to report status.
     *
     * @param value Percentage of completion between status reports.
     */
    void setUpdatePercent(const double value);

    /** Get how often to report status.
     *
     * @preturns Percentage of completion between status reports.
     */
    double getUpdatePercent(void) const;

    /** Set filename for output.
     *
     * @param filename Name of output file.
     */
    void setFilename(const char* filename);

    /** Get filename for output.
     *
     * @preturns Name of output file.
     */
    const char* getFilename(void) const;

    /// Open progress monitor.
    void open(void);

    /// Close progress monitor.
    void close(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    virtual
    void _open(void) = 0;

    /// Close progress monitor.
    virtual
    void _close(void) = 0;

    /** Compute finish time.
     * @param[in] percentComplete Percent of simulation complete.
     * @param[in] now Current data/time.
     * @param[in] startTime Date/time at start of simulation.
     * @returns String with estimated finish time.
     */
    std::string _calcFinishTime(const double percentComplete,
                                const time_t& now,
                                const time_t& startTime);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    double _updatePercent; ///< Percentage of completion between status reports.
    std::string _filename; ///< Name of output file.
    time_t _startTime;
    long _iUpdate; /// Current update step.
    bool _isMaster; ///< Is master process.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ProgressMonitor(const ProgressMonitor&); ///< Not implemented.
    const ProgressMonitor& operator=(const ProgressMonitor&); ///< Not implemented.

}; // ProgressMonitor

#endif // pylith_problems_progressmonitor_hh

// End of file
