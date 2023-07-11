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

#include "TestDataWriterHDF5ExtSubmesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtSubmesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5ExtSubmesh_Data* Tri(void);

    static TestDataWriterHDF5ExtSubmesh_Data* Quad(void);

    static TestDataWriterHDF5ExtSubmesh_Data* Tet(void);

    static TestDataWriterHDF5ExtSubmesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tri::testOpenClose", "[DataWriter][HDF5Ext][Submesh][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tri::testWriteVertexField", "[DataWriter][HDF5Ext][Submesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tri::testWriteCellField", "[DataWriter][HDF5Ext][Submesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtSubmesh::Quad::testOpenClose", "[DataWriter][HDF5Ext][Submesh][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Quad::testWriteVertexField", "[DataWriter][HDF5Ext][Submesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Quad::testWriteCellField", "[DataWriter][HDF5Ext][Submesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tet::testOpenClose", "[DataWriter][HDF5Ext][Submesh][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tet::testWriteVertexField", "[DataWriter][HDF5Ext][Submesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Tet::testWriteCellField", "[DataWriter][HDF5Ext][Submesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtSubmesh::Hex::testOpenClose", "[DataWriter][HDF5Ext][Submesh][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Hex::testWriteVertexField", "[DataWriter][HDF5Ext][Submesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtSubmesh::Hex::testWriteCellField", "[DataWriter][HDF5Ext][Submesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtSubmesh(pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtSubmesh_Data* data = new TestDataWriterHDF5ExtSubmesh_Data();assert(data);

    data->opencloseFilename = "tri3_surf.h5";
    data->vertexFilename = "tri3_surf_vertex.h5";
    data->cellFilename = "tri3_surf_cell.h5";

    TestDataWriterSubmesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtSubmesh_Data* data = new TestDataWriterHDF5ExtSubmesh_Data();assert(data);

    data->opencloseFilename = "quad4_surf.h5";
    data->vertexFilename = "quad4_surf_vertex.h5";
    data->cellFilename = "quad4_surf_cell.h5";

    TestDataWriterSubmesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtSubmesh_Data* data = new TestDataWriterHDF5ExtSubmesh_Data();assert(data);

    data->opencloseFilename = "tet4_surf.h5";
    data->vertexFilename = "tet4_surf_vertex.h5";
    data->cellFilename = "tet4_surf_cell.h5";

    TestDataWriterSubmesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Data*
pylith::meshio::TestDataWriterHDF5ExtSubmesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtSubmesh_Data* data = new TestDataWriterHDF5ExtSubmesh_Data();assert(data);

    data->opencloseFilename = "hex8_surf.h5";
    data->vertexFilename = "hex8_surf_vertex.h5";
    data->cellFilename = "hex8_surf_cell.h5";

    TestDataWriterSubmesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
