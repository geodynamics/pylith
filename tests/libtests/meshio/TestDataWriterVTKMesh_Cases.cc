// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKMesh_Cases {
public:

    // Data factory methods
    static TestDataWriterVTKMesh_Data* Tri(void);

    static TestDataWriterVTKMesh_Data* Quad(void);

    static TestDataWriterVTKMesh_Data* Tet(void);

    static TestDataWriterVTKMesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterVTKMesh::testAccessors", "[DataWriter][VTK][testAccessors]") {
    pylith::meshio::TestDataWriterVTKMesh::testAccessors();
}
TEST_CASE("TestDataWriterVTKMesh::testVtkFilename", "[DataWriter][VTK][testVtkFilename]") {
    pylith::meshio::TestDataWriterVTKMesh::testVtkFilename();
}

TEST_CASE("TestDataWriterVTKMesh::Tri::testTimeStep", "[DataWriter][VTK][Mesh][Tri][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tri()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMesh::Tri::testWriteVertexField", "[DataWriter][VTK][Mesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMesh::Tri::testWriteCellField", "[DataWriter][VTK][Mesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMesh::Quad::testTimeStep", "[DataWriter][VTK][Mesh][Quad][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Quad()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMesh::Quad::testWriteVertexField", "[DataWriter][VTK][Mesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMesh::Quad::testWriteCellField", "[DataWriter][VTK][Mesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMesh::Tet::testTimeStep", "[DataWriter][VTK][Mesh][Tet][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tet()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMesh::Tet::testWriteVertexField", "[DataWriter][VTK][Mesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMesh::Tet::testWriteCellField", "[DataWriter][VTK][Mesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMesh::Hex::testTimeStep", "[DataWriter][VTK][Mesh][Hex][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Hex()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMesh::Hex::testWriteVertexField", "[DataWriter][VTK][Mesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMesh::Hex::testWriteCellField", "[DataWriter][VTK][Mesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMesh(pylith::meshio::TestDataWriterVTKMesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMesh_Data*
pylith::meshio::TestDataWriterVTKMesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMesh_Data* data = new TestDataWriterVTKMesh_Data();assert(data);

    data->timestepFilename = "tri3.vtk";
    data->vertexFilename = "tri3_vertex.vtk";
    data->cellFilename = "tri3_cell.vtk";

    TestDataWriterMesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMesh_Data*
pylith::meshio::TestDataWriterVTKMesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMesh_Data* data = new TestDataWriterVTKMesh_Data();assert(data);

    data->timestepFilename = "quad4.vtk";
    data->vertexFilename = "quad4_vertex.vtk";
    data->cellFilename = "quad4_cell.vtk";

    TestDataWriterMesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMesh_Data*
pylith::meshio::TestDataWriterVTKMesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMesh_Data* data = new TestDataWriterVTKMesh_Data();assert(data);

    data->timestepFilename = "tet4.vtk";
    data->vertexFilename = "tet4_vertex.vtk";
    data->cellFilename = "tet4_cell.vtk";

    TestDataWriterMesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMesh_Data*
pylith::meshio::TestDataWriterVTKMesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMesh_Data* data = new TestDataWriterVTKMesh_Data();assert(data);

    data->timestepFilename = "hex8.vtk";
    data->vertexFilename = "hex8_vertex.vtk";
    data->cellFilename = "hex8_cell.vtk";

    TestDataWriterMesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
