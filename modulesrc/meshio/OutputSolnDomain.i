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
 * @file modulesrc/meshio/OutputSolnDomain.i
 *
 * @brief Python interface to C++ OutputSolnDomain object.
 */

namespace pylith {
    namespace meshio {
        class OutputSolnDomain : public pylith::meshio::OutputSoln {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor.
            OutputSolnDomain(void);

            /// Destructor
            virtual ~OutputSolnDomain(void);

            // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////
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

        }; // OutputSolnDomain

    } // meshio
} // pylith

// End of file
