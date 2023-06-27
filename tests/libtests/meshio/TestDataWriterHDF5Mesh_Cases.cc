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

#include "TestDataWriterHDF5Mesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Mesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Mesh_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5Mesh_Data* Tri(void);

    static TestDataWriterHDF5Mesh_Data* Quad(void);

    static TestDataWriterHDF5Mesh_Data* Tet(void);

    static TestDataWriterHDF5Mesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5Mesh::testAccessors", "[DataWriter][HDF5][testAccessors]") {
    pylith::meshio::TestDataWriterHDF5Mesh::testAccessors();
}
TEST_CASE("TestDataWriterHDF5Mesh::testHDF5Filename", "[DataWriter][HDF5][testHdf5Filename]") {
    pylith::meshio::TestDataWriterHDF5Mesh::testHdf5Filename();
}

TEST_CASE("TestDataWriterHDF5Mesh::Tri::testOpenClose", "[DataWriter][HDF5][Mesh][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Mesh::Tri::testWriteVertexField", "[DataWriter][HDF5][Mesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Mesh::Tri::testWriteCellField", "[DataWriter][HDF5][Mesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Mesh::Quad::testOpenClose", "[DataWriter][HDF5][Mesh][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Mesh::Quad::testWriteVertexField", "[DataWriter][HDF5][Mesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Mesh::Quad::testWriteCellField", "[DataWriter][HDF5][Mesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Mesh::Tet::testOpenClose", "[DataWriter][HDF5][Mesh][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Mesh::Tet::testWriteVertexField", "[DataWriter][HDF5][Mesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Mesh::Tet::testWriteCellField", "[DataWriter][HDF5][Mesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Mesh::Hex::testOpenClose", "[DataWriter][HDF5][Mesh][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Mesh::Hex::testWriteVertexField", "[DataWriter][HDF5][Mesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Mesh::Hex::testWriteCellField", "[DataWriter][HDF5][Mesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Mesh(pylith::meshio::TestDataWriterHDF5Mesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Mesh_Data*
pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Mesh_Data* data = new TestDataWriterHDF5Mesh_Data();assert(data);

    data->opencloseFilename = "tri3.h5";
    data->vertexFilename = "tri3_vertex.h5";
    data->cellFilename = "tri3_cell.h5";

    TestDataWriterMesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Mesh_Data*
pylith::meshio::TestDataWriterHDF5Mesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Mesh_Data* data = new TestDataWriterHDF5Mesh_Data();assert(data);

    data->opencloseFilename = "quad4.h5";
    data->vertexFilename = "quad4_vertex.h5";
    data->cellFilename = "quad4_cell.h5";

    TestDataWriterMesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Mesh_Data*
pylith::meshio::TestDataWriterHDF5Mesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Mesh_Data* data = new TestDataWriterHDF5Mesh_Data();assert(data);

    data->opencloseFilename = "tet4.h5";
    data->vertexFilename = "tet4_vertex.h5";
    data->cellFilename = "tet4_cell.h5";

    TestDataWriterMesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Mesh_Data*
pylith::meshio::TestDataWriterHDF5Mesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Mesh_Data* data = new TestDataWriterHDF5Mesh_Data();assert(data);

    data->opencloseFilename = "hex8.h5";
    data->vertexFilename = "hex8_vertex.h5";
    data->cellFilename = "hex8_cell.h5";

    TestDataWriterMesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
