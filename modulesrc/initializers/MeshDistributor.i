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
        class MeshDistributor : public pylith::initializers::InitializePhase {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            MeshDistributor(void);

            /// Default destructor
            ~MeshDistributor(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set distributor.
             *
             * @param[in] distributor Mesh distributor;
             */
            void setDistributor(pylith::topology::Distributor* distributor);

            /** Run initialization phase.
             *
             * @param[in] mesh Input mesh.
             * @param[in] problem Problem specification.
             * @returns Mesh after initialization phase.
             */
            pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                        const pylith::problems::Problem& problem);

        }; // MeshDistributor

    } // initializers
} // pylith

// End of file
