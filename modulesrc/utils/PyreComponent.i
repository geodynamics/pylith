// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
