// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


/**
 * @file modulesrc/topology/RefineUniform.hh
 *
 * @brief Python interface to C++ PyLith RefineUniform object.
 */

namespace pylith {
    namespace topology {
        class pylith::topology::RefineMesh : public pylith::utils::PyreComponent {
            // PUBLIC MEMBERS ///////////////////////////////////////////////////////
          public:
            /// Constructor
            RefineMesh(void);

            /// Destructor
            virtual ~RefineMesh(void);

            /// Deallocate data structures.
            void deallocate(void);

            /** Refine mesh.
             *
             * @param mesh Mesh to refine.
             * @returns Mesh after refinement.
             */
            virtual pylith::topology::Mesh* refine(const pylith::topology::Mesh& mesh) = 0;

        }; // RefineMesh

    } // namespace topology
} // namespace pylith


// End of file
