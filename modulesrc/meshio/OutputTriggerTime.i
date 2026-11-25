// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputTriggerTime : public pylith::meshio::OutputTrigger {
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

    } // namespace meshio
} // namespace pylith


// End of file
