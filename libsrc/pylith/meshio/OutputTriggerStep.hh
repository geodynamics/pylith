// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputTriggerStep.hh
 *
 * @brief Base decision on whether to write output based on number of time steps since more recent write.
 */

#if !defined(pylith_meshio_outputtriggerstep_hh)
#define pylith_meshio_outputtriggerstep_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/meshio/OutputTrigger.hh" // ISA OutputTrigger

// OutputTrigger --------------------------------------------------------
/// Abstract base class for output trigger.
class pylith::meshio::OutputTriggerStep : public pylith::meshio::OutputTrigger {
    friend class TestOutputTriggerStep;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputTriggerStep(void);

    /// Destructor
    virtual ~OutputTriggerStep(void);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex);

    /** Set number of time steps to skip between writes.
     *
     * @param[in] Number of time steps to skip between writes.
     */
    void numTimeStepsSkip(const int value);

    /** Get number of time steps to skip between writes.
     *
     * @returns Number of time steps to skip between writes.
     */
    int numTimeStepsSkip(void) const;


    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    PylithInt _numTimeStepsSkip; ///< Number of time steps to skip between writes.
    PylithInt _timeStepWrote; ///< Time step when data was previously written.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputTriggerStep(const OutputTriggerStep&);   ///< Not implemented.
    const OutputTriggerStep& operator=(const OutputTriggerStep&);   ///< Not implemented

}; // OutputTriggerStep

#endif // pylith_meshio_outputtriggerstep_hh


// End of file
