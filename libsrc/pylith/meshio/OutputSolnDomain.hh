// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include "pylith/meshio/OutputSoln.hh" // ISA OutputSoln
#include "pylith/problems/problemsfwd.hh" // HASA Problem

class pylith::meshio::OutputSolnDomain : public pylith::meshio::OutputSoln {
    friend class TestOutputSolnDomain; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    OutputSolnDomain(void);

    /// Destructor
    virtual ~OutputSolnDomain(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void _writeSolnStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputSolnDomain(const OutputSolnDomain&); ///< Not implemented.
    const OutputSolnDomain& operator=(const OutputSolnDomain&); ///< Not implemented

}; // OutputSolnDomain

// End of file
