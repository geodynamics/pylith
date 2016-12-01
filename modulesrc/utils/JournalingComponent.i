// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/utils/JournalingComponent.i
 *
 * @brief Python interface to C++ JournalingComponent.
 */


namespace pylith {
    namespace utils {

        class JournalingComponent
        { // JournalingComponent

        // PUBLIC MEMBERS /////////////////////////////////////////////////
public:

        /// Constructor
        JournalingComponent(void);

        /// Destructor
        ~JournalingComponent(void);

        /** Set name of journal.
         *
         * @param value Name of journal.
         */
        void name(const char* value);

        /** Get name of journal.
         *
         * @returns Name of journal.
         */
        const char* name(void) const;

        /// Setup journaling..
        void initialize(void);

        /** Get debug journal.
         *
         * @returns Debugging journal.
         */
        journal::debug_t& debug(void);

        /** Get info journal.
         *
         * @returns Info journal.
         */
        journal::info_t& info(void);

        /** Get warning journal.
         *
         * @returns Warning journal.
         */
        journal::warning_t& warning(void);

        /** Get error journal.
         *
         * @returns Error journal.
         */
        journal::error_t& error(void);

        }; // JournalingComponent

    } // utils
} // pylith


// End of file
