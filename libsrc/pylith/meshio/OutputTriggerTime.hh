// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file libsrc/meshio/OutputTriggerTime.hh
 *
 * @brief Base decision on whether to write output based on time since more recent write.
 */

#if !defined(pylith_meshio_outputtriggertime_hh)
#define pylith_meshio_outputtriggertime_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/meshio/OutputTrigger.hh" // ISA OutputTrigger

class pylith::meshio::OutputTriggerTime : public pylith::meshio::OutputTrigger {
    friend class TestOutputTriggerTime; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    OutputTriggerTime(void);

    /// Destructor
    virtual ~OutputTriggerTime(void);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex);

    /** Set elapsed time between writes.
     *
     * @param[in] Elapsed time between writes.
     */
    void setTimeSkip(const double value);

    /** Get elapsed time between writes.
     *
     * @returns Elapsed time between writes.
     */
    double getTimeSkip(void) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PylithReal _timeSkip; ///< Elapsed (dimensional) time between writes.
    PylithReal _timeNondimWrote; ///< Time (nondimensional) when data was previously writtern.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputTriggerTime(const OutputTriggerTime&); ///< Not implemented.
    const OutputTriggerTime& operator=(const OutputTriggerTime&); ///< Not implemented

};

// OutputTriggerTime

#endif // pylith_meshio_outputtriggertime_hh

// End of file
