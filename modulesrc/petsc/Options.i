// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/petsc/Options.i
 *
 * @brief Python interface to C++ Defaults.
 */

namespace pylith {
    namespace petsc {
        class Defaults: public pylith::utils::GenericComponent {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const int NONE;
            static const int MONITORS;
            static const int SOLVER;
            static const int PARALLEL;
            static const int INITIAL_GUESS;
            static const int TESTING;
            static const int COLLECTIVE_IO;
            static const int TS_ADAPTIVE;
            static const int TS_IMPULSE;

            // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////
private:

            Defaults(void); ///< Not implemented
            Defaults(const Defaults &); ///< Not implemented.
            const Defaults& operator = (const Defaults&); ///< Not implemented

        }; // class Defaults

    } // utils
} // pylith

// End of file
