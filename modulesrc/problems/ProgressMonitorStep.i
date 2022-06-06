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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/problems/ProgressMonitorStep.i
 *
 * Python interface to C++ class ProgressMonitorStep.
 */

namespace pylith {
    namespace problems {
        class ProgressMonitorStep: public pylith::problems::ProgressMonitor {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            ProgressMonitorStep(void);

            /// Destructor
            virtual ~ProgressMonitorStep(void);

            /** Update progress.
             *
             * @param[in] current Current step.
             * @param[in] start Starting step.
             * @param[in] stop Ending step.
             */
            void update(const size_t current,
                        const size_t start,
                        const size_t stop);

            // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////
protected:

            /// Open progress monitor.
            void _open(void);

            /// Close progress monitor.
            void _close(void);

            /** Update progress.
             *
             * @param[in] setp Current step.
             * @param[in] now Current date/time.
             * @param[in] percentComplete Percent completed
             * @param[in] finished Time stamp of estimated finish.
             */
            void _update(const size_t step,
                         const time_t& now,
                         const double percentComplete,
                         const char* finished);

        }; // class ProgressMonitorStep

    } // problems
} // pylith

// End of file
