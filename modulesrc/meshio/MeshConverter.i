// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace meshio {
        class MeshConverter : public pylith::utils::GenericComponent {
            // PUBLIC METHODS /////////////////////////////////////////////////
          public:
            /** Convert mesh format.
             *
             * @param[in] writer Mesh writer.
             * @param[in] reader Mesh reader.
             * @param[in] checkTopology Check topology of input mesh.
             */
            static void convert(pylith::meshio::MeshIO* writer,
                                pylith::meshio::MeshIO* reader,
                                const bool checkTopology = true);

        }; // MeshConverter

    } // namespace meshio
} // namespace pylith


// End of file
