// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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

    /** Update progress.
     *
     * @param[in] current Current time.
     * @param[in] start Starting time.
     * @param[in] stop Ending time.
     */
    void update(const double current,
                const double start,
                const double stop);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Open progress monitor.
    void _open(void);

    /// Close progress monitor.
    void _close(void);

    /** Update progress.
     *
     * @param[in] t Current time.
     * @param[in] now Current date/time.
     * @param[in] percentComplete Percent completed
     * @param[in] finished Time stamp of estimated finish.
     */
    void _update(const double t,
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
