// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/meshio/meshiofwd.hh" // USES MeshIO

#include "pylith/meshio/MeshBuilder.hh" // HOLDSA Topology, Geometry

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

    pylith::meshio::MeshBuilder::Topology* topology;
    pylith::meshio::MeshBuilder::Geometry* geometry;

    PylithInt* materialIds; ///< Pointer to cell material identifiers

    PylithInt* vertexGroups; ///< Array of pointers to indices of points in vertex groups
    PylithInt* vertexGroupSizes; ///< Array of sizes of each vertex group
    PylithInt* vertexGroupTags; ///< Array of label values (tags) for each vertex group.
    char** vertexGroupNames; ///< Array of vertex group names
    size_t numVertexGroups; ///< Number of vertex groups

    PylithInt* faceGroups; ///< Array of pointers to indices of points in face groups
    PylithInt* faceGroupSizes; ///< Array of sizes of each face group
    PylithInt* faceGroupTags; ///< Array of label values (tags) for each face group.
    char** faceGroupNames; ///< Array of face group names
    size_t numFaceGroups; ///< Number of face groups
    size_t numFaceVertices; ///< Number of vertices on a cell face.

    bool useIndexZero; ///< Indices start with 0 if true, 1 if false

    std::string filename;
};

// End of file
