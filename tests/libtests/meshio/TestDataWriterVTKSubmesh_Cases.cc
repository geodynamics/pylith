// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestDataWriterVTKSubmesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterVTKSubmesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKSubmesh_Cases {
public:

    // Data factory methods
    static TestDataWriterVTKSubmesh_Data* Tri(void);

    static TestDataWriterVTKSubmesh_Data* Quad(void);

    static TestDataWriterVTKSubmesh_Data* Tet(void);

    static TestDataWriterVTKSubmesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterVTKMesh::Tri::testTimeStep", "[DataWriter][VTK][Submesh][Tri][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tri()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKSubmesh::Tri::testWriteVertexField", "[DataWriter][VTK][Submesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKSubmesh::Tri::testWriteCellField", "[DataWriter][VTK][Submesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKSubmesh::Quad::testTimeStep", "[DataWriter][VTK][Submesh][Quad][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Quad()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKSubmesh::Quad::testWriteVertexField", "[DataWriter][VTK][Submesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKSubmesh::Quad::testWriteCellField", "[DataWriter][VTK][Submesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKSubmesh::Tet::testTimeStep", "[DataWriter][VTK][Submesh][Tet][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tet()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKSubmesh::Tet::testWriteVertexField", "[DataWriter][VTK][Submesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKSubmesh::Tet::testWriteCellField", "[DataWriter][VTK][Submesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKSubmesh::Hex::testTimeStep", "[DataWriter][VTK][Submesh][Hex][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Hex()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKSubmesh::Hex::testWriteVertexField", "[DataWriter][VTK][Submesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKSubmesh::Hex::testWriteCellField", "[DataWriter][VTK][Submesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKSubmesh(pylith::meshio::TestDataWriterVTKSubmesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKSubmesh_Data*
pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKSubmesh_Data* data = new TestDataWriterVTKSubmesh_Data();assert(data);

    data->timestepFilename = "tri3_surf.vtu";
    data->vertexFilename = "tri3_surf_vertex.vtu";
    data->cellFilename = "tri3_surf_cell.vtu";

    TestDataWriterSubmesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKSubmesh_Data*
pylith::meshio::TestDataWriterVTKSubmesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKSubmesh_Data* data = new TestDataWriterVTKSubmesh_Data();assert(data);

    data->timestepFilename = "quad4_surf.vtu";
    data->vertexFilename = "quad4_surf_vertex.vtu";
    data->cellFilename = "quad4_surf_cell.vtu";

    TestDataWriterSubmesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKSubmesh_Data*
pylith::meshio::TestDataWriterVTKSubmesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKSubmesh_Data* data = new TestDataWriterVTKSubmesh_Data();assert(data);

    data->timestepFilename = "tet4_surf.vtu";
    data->vertexFilename = "tet4_surf_vertex.vtu";
    data->cellFilename = "tet4_surf_cell.vtu";

    TestDataWriterSubmesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKSubmesh_Data*
pylith::meshio::TestDataWriterVTKSubmesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKSubmesh_Data* data = new TestDataWriterVTKSubmesh_Data();assert(data);

    data->timestepFilename = "hex8_surf.vtu";
    data->vertexFilename = "hex8_surf_vertex.vtu";
    data->cellFilename = "hex8_surf_cell.vtu";

    TestDataWriterSubmesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
