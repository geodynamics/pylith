// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/meshio/MeshConverter.i
 *
 * @brief Python interface to C++ MeshConverter object.
 */

namespace pylith {
    namespace meshio {
        class MeshConverter: public pylith::utils::GenericComponent {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /** Convert mesh format.
             *
             * @param[in] writer Mesh writer.
             * @param[in] reader Mesh reader.
             * @param[in] checkTopology Check topology of input mesh.
             */
            static
            void convert(pylith::meshio::MeshIO* writer,
                         pylith::meshio::MeshIO* reader,
                         const bool checkTopology=true);

        }; // MeshConverter

    } // meshio
} // pylith

// End of file
