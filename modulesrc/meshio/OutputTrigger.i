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
            virtual bool shouldWrite(const PylithReal t,
                                     const PylithInt tindex) = 0;

        }; // OutputTrigger

    } // namespace meshio
} // namespace pylith


// End of file
