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
 * @file libsrc/meshio/OutputTrigger.hh
 *
 * @brief Abstract base class for output trigger.
 */

#if !defined(pylith_meshio_outputtrigger_hh)
#define pylith_meshio_outputtrigger_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh"

#include "pylith/utils/types.hh" // USE PylithInt, PylithReal

// OutputTrigger --------------------------------------------------------
/// Abstract base class for output trigger.
class pylith::meshio::OutputTrigger : public pylith::utils::PyreComponent {
    friend class TestOutputTrigger;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputTrigger(void);

    /// Destructor
    virtual ~OutputTrigger(void);

    /** Check whether we want to write output at time t.
     *
     * @param[in] t Time of proposed write.
     * @param[in] tindex Inxex of current time step.
     * @returns True if output should be written at time t, false otherwise.
     */
    virtual
    bool shouldWrite(const PylithReal t,
                     const PylithInt tindex) = 0;


    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputTrigger(const OutputTrigger&);   ///< Not implemented.
    const OutputTrigger& operator=(const OutputTrigger&);   ///< Not implemented

}; // OutputTrigger

#endif // pylith_meshio_outputtrigger_hh


// End of file
