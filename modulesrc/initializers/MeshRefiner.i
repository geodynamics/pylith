// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


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

    } // namespace initializers
} // namespace pylith


// End of file
