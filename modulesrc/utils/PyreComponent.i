// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file modulesrc/utils/PyreComponent.i
 *
 * @brief Python interface to C++ PyreComponent.
 */


namespace pylith {
    namespace utils {

        class PyreComponent {

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
        void setName(const char* value);

        /** Get name of component.
         *
         * @returns Name of component.
         */
        const char* getName(void) const;

        /** Set component identifier (identifies object in component hierarchy).
         *
         * @param value Component identifier.
         */
        void setIdentifier(const char* value);

        /** Get component identifier (identifies object in component hierarchy).
         *
         * @returns Component identifier.
         */
        const char* getIdentifier(void) const;

        }; // PyreComponent

    } // utils
} // pylith


// End of file
