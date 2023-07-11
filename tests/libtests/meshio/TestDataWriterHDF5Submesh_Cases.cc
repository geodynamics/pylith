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

#include "TestDataWriterHDF5Submesh.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Submesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Submesh_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5Submesh_Data* Tri(void);

    static TestDataWriterHDF5Submesh_Data* Quad(void);

    static TestDataWriterHDF5Submesh_Data* Tet(void);

    static TestDataWriterHDF5Submesh_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5Submesh::Tri::testOpenClose", "[DataWriter][HDF5][Submesh][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Submesh::Tri::testWriteVertexField", "[DataWriter][HDF5][Submesh][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Submesh::Tri::testWriteCellField", "[DataWriter][HDF5][Submesh][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Submesh::Quad::testOpenClose", "[DataWriter][HDF5][Submesh][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Submesh::Quad::testWriteVertexField", "[DataWriter][HDF5][Submesh][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Submesh::Quad::testWriteCellField", "[DataWriter][HDF5][Submesh][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Submesh::Tet::testOpenClose", "[DataWriter][HDF5][Submesh][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Submesh::Tet::testWriteVertexField", "[DataWriter][HDF5][Submesh][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Submesh::Tet::testWriteCellField", "[DataWriter][HDF5][Submesh][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Submesh::Hex::testOpenClose", "[DataWriter][HDF5][Submesh][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Submesh::Hex::testWriteVertexField", "[DataWriter][HDF5][Submesh][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Submesh::Hex::testWriteCellField", "[DataWriter][HDF5][Submesh][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Submesh(pylith::meshio::TestDataWriterHDF5Submesh_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Submesh_Data*
pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Submesh_Data* data = new TestDataWriterHDF5Submesh_Data();assert(data);

    data->opencloseFilename = "tri3_surf.h5";
    data->vertexFilename = "tri3_surf_vertex.h5";
    data->cellFilename = "tri3_surf_cell.h5";

    TestDataWriterSubmesh::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Submesh_Data*
pylith::meshio::TestDataWriterHDF5Submesh_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Submesh_Data* data = new TestDataWriterHDF5Submesh_Data();assert(data);

    data->opencloseFilename = "quad4_surf.h5";
    data->vertexFilename = "quad4_surf_vertex.h5";
    data->cellFilename = "quad4_surf_cell.h5";

    TestDataWriterSubmesh::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Submesh_Data*
pylith::meshio::TestDataWriterHDF5Submesh_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Submesh_Data* data = new TestDataWriterHDF5Submesh_Data();assert(data);

    data->opencloseFilename = "tet4_surf.h5";
    data->vertexFilename = "tet4_surf_vertex.h5";
    data->cellFilename = "tet4_surf_cell.h5";

    TestDataWriterSubmesh::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Submesh_Data*
pylith::meshio::TestDataWriterHDF5Submesh_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Submesh_Data* data = new TestDataWriterHDF5Submesh_Data();assert(data);

    data->opencloseFilename = "hex8_surf.h5";
    data->vertexFilename = "hex8_surf_vertex.h5";
    data->cellFilename = "hex8_surf_cell.h5";

    TestDataWriterSubmesh::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
