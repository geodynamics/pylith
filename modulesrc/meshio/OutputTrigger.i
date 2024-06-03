// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/meshio/OutputTrigger.i
 *
 * @brief Python interface to C++ OutputTrigger object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputTrigger: public pylith::utils::PyreComponent {
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
