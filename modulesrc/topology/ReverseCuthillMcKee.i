// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


/**
 * @file modulesrc/topology/ReverseCuthillMcKee.hh
 *
 * @brief Python interface to C++ PyLith ReverseCuthillMcKee object.
 */

namespace pylith {
    namespace topology {
        // ReverseCuthillMcKee ----------------------------------------------
        class ReverseCuthillMcKee { // ReverseCuthillMcKee
                                    // PUBLIC METHODS /////////////////////////////////////////////////
          public:
            /** Reorder vertices and cells of mesh using PETSc routines
             * implementing reverse Cuthill-McKee algorithm.
             *
             * @param mesh PyLith finite-element mesh.
             */
            static pylith::topology::Mesh* reorder(const pylith::topology::Mesh& mesh);

        }; // ReverseCuthillMcKee

    } // namespace topology
} // namespace pylith


// End of file
