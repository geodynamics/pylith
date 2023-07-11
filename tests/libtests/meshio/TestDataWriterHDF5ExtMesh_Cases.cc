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

#include "TestDataWriterHDF5ExtMesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtMesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtMesh_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5ExtMesh_Data* Tri(void);

    static TestDataWriterHDF5ExtMesh_Data* Quad(void);

    static TestDataWriterHDF5ExtMesh_Data* Tet(void);

    static TestDataWriterHDF5ExtMesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5ExtMesh::testAccessors", "[DataWriter][HDF5Ext][testAccessors]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh::testAccessors();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::testHDF5Filename", "[DataWriter][HDF5Ext][testHdf5Filename]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh::testHdf5Filename();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::testDatasetFilename", "[DataWriter][HDF5Ext][testDatasetFilename]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh::testDatasetFilename();
}

TEST_CASE("TestDataWriterHDF5ExtMesh::Tri::testOpenClose", "[DataWriter][HDF5Ext][Mesh][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Tri::testWriteVertexField", "[DataWriter][HDF5Ext][Mesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Tri::testWriteCellField", "[DataWriter][HDF5Ext][Mesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMesh::Quad::testOpenClose", "[DataWriter][HDF5Ext][Mesh][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Quad::testWriteVertexField", "[DataWriter][HDF5Ext][Mesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Quad::testWriteCellField", "[DataWriter][HDF5Ext][Mesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMesh::Tet::testOpenClose", "[DataWriter][HDF5Ext][Mesh][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Tet::testWriteVertexField", "[DataWriter][HDF5Ext][Mesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Tet::testWriteCellField", "[DataWriter][HDF5Ext][Mesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMesh::Hex::testOpenClose", "[DataWriter][HDF5Ext][Mesh][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Hex::testWriteVertexField", "[DataWriter][HDF5Ext][Mesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMesh::Hex::testWriteCellField", "[DataWriter][HDF5Ext][Mesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMesh(pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMesh_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMesh_Data* data = new TestDataWriterHDF5ExtMesh_Data();assert(data);

    data->opencloseFilename = "tri3.h5";
    data->vertexFilename = "tri3_vertex.h5";
    data->cellFilename = "tri3_cell.h5";

    TestDataWriterMesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMesh_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMesh_Data* data = new TestDataWriterHDF5ExtMesh_Data();assert(data);

    data->opencloseFilename = "quad4.h5";
    data->vertexFilename = "quad4_vertex.h5";
    data->cellFilename = "quad4_cell.h5";

    TestDataWriterMesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMesh_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMesh_Data* data = new TestDataWriterHDF5ExtMesh_Data();assert(data);

    data->opencloseFilename = "tet4.h5";
    data->vertexFilename = "tet4_vertex.h5";
    data->cellFilename = "tet4_cell.h5";

    TestDataWriterMesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMesh_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMesh_Data* data = new TestDataWriterHDF5ExtMesh_Data();assert(data);

    data->opencloseFilename = "hex8.h5";
    data->vertexFilename = "hex8_vertex.h5";
    data->cellFilename = "hex8_cell.h5";

    TestDataWriterMesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
