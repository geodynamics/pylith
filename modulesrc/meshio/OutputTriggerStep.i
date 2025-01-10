// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/meshio/OutputTriggerStep.i
 *
 * @brief Python interface to C++ OutputTriggerStep object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputTriggerStep: public pylith::meshio::OutputTrigger {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor
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

        }; // OutputTriggerStep

    } // meshio
} // pylith

// End of file
