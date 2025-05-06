// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

// ------------------------------------------------------------------------------------------------
class pylith::meshio::MeshConverter : public pylith::utils::GenericComponent {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
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

// End of file
