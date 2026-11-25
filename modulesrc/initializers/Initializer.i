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

    } // namespace initializers
} // namespace pylith


// End of file
