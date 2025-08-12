// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 4, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/topology/RefineUniform.hh
 *
 * @brief Python interface to C++ PyLith RefineUniform object.
 */

namespace pylith {
    namespace topology {
        class RefineUniform : public pylith::topology::RefineMesh {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            RefineUniform(void);

            /// Destructor
            ~RefineUniform(void);

            /// Deallocate data structures.
            void deallocate(void);

            /** Set number of levels of refinement.
             *
             * @param[in] numLevels Number of levels.
             */
            void setNumLevels(const size_t numLevels);

            /** Get number of levels of refinement.
             *
             * @returns Number of levels.
             */
            size_t getNumLevels(void) const;

            /** Refine mesh.
             *
             * @param mesh Mesh to refine.
             * @returns Mesh after refinement.
             */
            pylith::topology::Mesh* refine(const pylith::topology::Mesh& mesh);

        }; // RefineUniform

    } // topology
} // pylith

// End of file
