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
 * @file modulesrc/utils/PyreComponent.i
 *
 * @brief Python interface to C++ PyreComponent.
 */


namespace pylith {
    namespace utils {

        class PyreComponent
        { // PyreComponent

        // PUBLIC MEMBERS /////////////////////////////////////////////////
public:

        /// Constructor
        PyreComponent(void);

        /// Destructor
        ~PyreComponent(void);

        /** Set name of component.
         *
         * @param value Name of component.
         */
        void name(const char* value);

        /** Get name of component.
         *
         * @returns Name of component.
         */
        const char* name(void) const;

        /** Set component identifier (identifies object in component hierarchy).
         *
         * @param value Component identifier.
         */
        void identifier(const char* value);

        /** Get component identifier (identifies object in component hierarchy).
         *
         * @returns Component identifier.
         */
        const char* identifier(void) const;

        }; // PyreComponent

    } // utils
} // pylith


// End of file
