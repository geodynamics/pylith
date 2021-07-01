// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file modulesrc/problems/ProgressMonitorTime.i
 *
 * Python interface to C++ abstract base class ProgressMonitorTime.
 */

namespace pylith {
    namespace problems {
        class ProgressMonitorTime : public pylith::problems::ProgressMonitor {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
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

            // PROTECTED MEMBERS
            // ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

            /// Open progress monitor.
            void _open(void);

            /// Close progress monitor.
            void _close(void);

            /** Update progress.
             *
             * @param[in current Current step.
             * @param[in] now Current time.
             * @param[in] percentComplete Percent completed
             * @param[in] finished Time stamp of estimated finish.
             */
            void _update(const double current,
                         const time_t& now,
                         const double percentComplete,
                         const char* finished);

        }; // class ProgressMonitorTime

    } // problems
} // pylith

// End of file
