// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

namespace pylith {
    namespace initializers {
        class MeshRefiner : public pylith::initializers::InitializePhase {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            MeshRefiner(void);

            /// Default destructor
            ~MeshRefiner(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set refiner.
             *
             * @param[in] refiner Mesh refiner;
             */
            void setRefiner(pylith::topology::RefineMesh* refiner);

            /** Run initialization phase.
             *
             * @param[in] mesh Input mesh.
             * @param[in] problem Problem specification.
             * @returns Mesh after initialization phase.
             */
            pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                        const pylith::problems::Problem& problem);

        }; // MeshRefiner

    } // initializers
} // pylith

// End of file
