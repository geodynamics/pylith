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
 * @file modulesrc/meshio/OutputTrigger.i
 *
 * @brief Python interface to C++ OutputTrigger object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputTrigger : public pylith::utils::PyreComponent {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
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

        }; // OutputTrigger

    } // meshio
} // pylith

// End of file
