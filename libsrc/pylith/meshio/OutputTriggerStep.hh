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
