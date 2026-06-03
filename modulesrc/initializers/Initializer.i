// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2026, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

namespace pylith {
    namespace initializers {
        class Initializer {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            Initializer(void);

            /// Default destructor
            ~Initializer(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set phases.
             *
             * @param[in] phases Initialization phases.
             */
            void setPhases(pylith::initializers::InitializePhase* phases[],
                           const size_t numPhases);

            /** Run initialization phase.
             *
             * @param[in] problem Problem specification.
             * @returns Mesh after initialization phase.
             */
            pylith::topology::Mesh* runPhases(const pylith::problems::Problem& problem);

        }; // Initializer

    } // initializers
} // pylith

// End of file
