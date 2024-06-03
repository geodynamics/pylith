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
 * @file modulesrc/meshio/OutputTriggerTime.i
 *
 * @brief Python interface to C++ OutputTriggerTime object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputTriggerTime: public pylith::meshio::OutputTrigger {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor
            OutputTriggerTime(void);

            /// Destructor
            ~OutputTriggerTime(void);

            /** Check whether we want to write output at time t.
             *
             * @param[in] t Time of proposed write.
             * @param[in] tindex Inxex of current time step.
             * @returns True if output should be written at time t, false otherwise.
             */
            bool shouldWrite(const PylithReal t,
                             const PylithInt tindex);

            /** Set elapsed time between writes.
             *
             * @param[in] Elapsed time between writes.
             */
            void setTimeSkip(const double value);

            /** Get elapsed time between writes.
             *
             * @returns Elapsed time between writes.
             */
            double getTimeSkip(void) const;

        }; // OutputTriggerTime

    } // meshio
} // pylith

// End of file
