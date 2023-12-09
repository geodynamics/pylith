// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/utils/PetscOptions.i
 *
 * @brief Python interface to C++ PetscDefaults.
 */

namespace pylith {
    namespace utils {
        class pylith::utils::PetscDefaults: public pylith::utils::GenericComponent {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const int NONE;
            static const int MONITORS;
            static const int SOLVER;
            static const int PARALLEL;
            static const int INITIAL_GUESS;
            static const int TESTING;

            // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////
private:

            PetscDefaults(void); ///< Not implemented
            PetscDefaults(const PetscDefaults &); ///< Not implemented.
            const PetscDefaults& operator = (const PetscDefaults&); ///< Not implemented

        }; // class PetscDefaults

    } // utils
} // pylith

// End of file
