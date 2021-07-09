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
 * @file libsrc/meshio/OutputTriggerStep.hh
 *
 * @brief Base decision on whether to write output based on number of steps since more recent write.
 */

#if !defined(pylith_meshio_outputtriggerstep_hh)
#define pylith_meshio_outputtriggerstep_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/meshio/OutputTrigger.hh" // ISA OutputTrigger

class pylith::meshio::OutputTriggerStep : public pylith::meshio::OutputTrigger {
    friend class TestOutputTriggerStep; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Construlibctor
    OutputTriggerStep(void);

    /// Destructor
    ~OutputTriggerStep(void);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex);

    /** Set number of steps to skip between writes.
     *
     * @param[in] Number of steps to skip between writes.
     */
    void setNumStepsSkip(const int value);

    /** Get number of steps to skip between writes.
     *
     * @returns Number of steps to skip between writes.
     */
    int getNumStepsSkip(void) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PylithInt _numStepsSkip; ///< Number of steps to skip between writes.
    PylithInt _stepWrote; ///< Step when data was previously written.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputTriggerStep(const OutputTriggerStep&); ///< Not implemented.
    const OutputTriggerStep& operator=(const OutputTriggerStep&); ///< Not implemented

};

// OutputTriggerStep

#endif // pylith_meshio_outputtriggerstep_hh

// End of file
