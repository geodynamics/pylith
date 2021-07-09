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
 * @file modulesrc/meshio/OutputSolnBoundary.i
 *
 * @brief Python interface to C++ OutputSolnBoundary object.
 */

namespace pylith {
    namespace meshio {
        class OutputSolnBoundary : public pylith::meshio::OutputSoln {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor.
            OutputSolnBoundary(void);

            /// Destructor
            virtual ~OutputSolnBoundary(void);

            /** Set label identifier for subdomain.
             *
             * @param[in] value Label of subdomain.
             */
            void setLabel(const char* value);

            /** Verify configuration.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

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

        }; // OutputSolnBoundary

    } // meshio
} // pylith

// End of file
