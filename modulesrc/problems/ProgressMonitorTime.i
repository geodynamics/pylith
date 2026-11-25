// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace problems {
        class ProgressMonitorTime : public pylith::problems::ProgressMonitor {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
          public:
            /// Constructor
            ProgressMonitorTime(void);

            /// Destructor
            virtual ~ProgressMonitorTime(void);

            /** Set unit for simulation time in output.
             *
             * @param[in] Unit of time.
             */
            void setTimeUnit(const char* value);

            /** Set unit for simulation time in output.
             *
             * @param[in] Unit of time.
             */
            const char* getTimeUnit(void) const;

            /** Update progress.
             *
             * @param[in] current Current time.
             * @param[in] start Starting time.
             * @param[in] stop Ending time.
             */
            void update(const double current,
                        const double start,
                        const double stop);

            // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////
          protected:
            /// Open progress monitor.
            void _open(void);

            /// Close progress monitor.
            void _close(void);

            /** Update progress.
             *
             * @param[in] t Current time.
             * @param[in] now Current date/time.
             * @param[in] percentComplete Percent completed
             * @param[in] finished Time stamp of estimated finish.
             */
            void _update(const double t,
                         const time_t& now,
                         const double percentComplete,
                         const char* finished);

        }; // class ProgressMonitorTime

    } // namespace problems
} // namespace pylith


// End of file
