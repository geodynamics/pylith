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
 * @file libsrc/meshio/OutputTrigger.hh
 *
 * @brief Abstract base class for output trigger.
 */

#if !defined(pylith_meshio_outputtrigger_hh)
#define pylith_meshio_outputtrigger_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh"

#include "pylith/utils/types.hh" // USE PylithInt, PylithReal

class pylith::meshio::OutputTrigger : public pylith::utils::PyreComponent {
    friend class TestOutputTrigger; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    OutputTrigger(void);

    /// Destructor
    virtual ~OutputTrigger(void);

    /** Set time scale for nondimensionalizing time.
     *
     * @param[in] value Time scale.
     */
    void setTimeScale(const PylithReal value);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    virtual
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex) = 0;

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _timeScale; ///< Time scale for nondimensionalizing time.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputTrigger(const OutputTrigger&); ///< Not implemented.
    const OutputTrigger& operator=(const OutputTrigger&); ///< Not implemented

};

// OutputTrigger

#endif // pylith_meshio_outputtrigger_hh

// End of file
