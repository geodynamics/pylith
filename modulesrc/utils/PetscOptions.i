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
