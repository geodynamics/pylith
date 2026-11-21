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
        class InitializePhase {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////
          public:
            /// Default constructor.
            InitializePhase(void);

            /// Default destructor
            ~InitializePhase(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Run initialization phase.
             *
             * @param[in] mesh Input mesh.
             * @param[in] problem Problem specification.
             * @returns Mesh after initialization phase.
             */
            virtual pylith::topology::Mesh* run(pylith::topology::Mesh* mesh,
                                                const pylith::problems::Problem& problem) = 0;

        }; // InitializePhase

    } // namespace initializers
} // namespace pylith


// End of file
