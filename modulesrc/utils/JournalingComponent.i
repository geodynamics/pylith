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

        }; // JournalingComponent

    } // utils
} // pylith


// End of file
