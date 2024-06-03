// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/meshio/meshiofwd.hh" // USES MeshIO

namespace pylith {
    namespace meshio {
        class TestMeshIO;
        class TestMeshIO_Data; // test data
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIO : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIO(TestMeshIO_Data* data);

    /// Tear down testing data.
    ~TestMeshIO(void);

    /// Get simple mesh for testing output.
    void _createMesh(void);

    /// Check values in mesh against data.
    void _checkVals(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestMeshIO_Data* _data; ///< Test data.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestMeshIO

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIO_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestMeshIO_Data(void);

    /// Destructor
    ~TestMeshIO_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    PylithInt numVertices; ///< Number of vertices
    PylithInt spaceDim; ///< Number of dimensions in vertex coordinates
    PylithInt numCells; ///< Number of cells
    PylithInt cellDim; ///< Number of dimensions associated with cell
    PylithInt numCorners; ///< Number of vertices in cell

    PylithScalar* vertices; ///< Pointer to coordinates of vertices
    PylithInt* cells; ///< Pointer to indices of vertices in cells
    PylithInt* materialIds; ///< Pointer to cell material identifiers

    PylithInt* groups; ///< Array of pointers to indices of points in groups
    PylithInt* groupSizes; ///< Array of sizes of each group
    PylithInt* groupTags; ///< Array of label values (tags) for each group.
    char** groupNames; ///< Array of group names
    char** groupTypes; ///< Array of group types
    PylithInt numGroups; ///< Number of groups

    bool useIndexZero; ///< Indices start with 0 if true, 1 if false

    std::string filename;
};

// End of file
