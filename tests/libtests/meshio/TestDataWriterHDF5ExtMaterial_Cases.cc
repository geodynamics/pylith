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

#include "TestDataWriterHDF5ExtMaterial.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtMaterial_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5ExtMaterial_Data* Tri(void);

    static TestDataWriterHDF5ExtMaterial_Data* Quad(void);

    static TestDataWriterHDF5ExtMaterial_Data* Tet(void);

    static TestDataWriterHDF5ExtMaterial_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5ExtMaterial::Tri::testOpenClose", "[DataWriter][HDF5Ext][Material][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Tri::testWriteVertexField", "[DataWriter][HDF5Ext][Material][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Tri::testWriteCellField", "[DataWriter][HDF5Ext][Material][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMaterial::Quad::testOpenClose", "[DataWriter][HDF5Ext][Material][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Quad::testWriteVertexField", "[DataWriter][HDF5Ext][Material][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Quad::testWriteCellField", "[DataWriter][HDF5Ext][Material][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMaterial::Tet::testOpenClose", "[DataWriter][HDF5Ext][Material][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Tet::testWriteVertexField", "[DataWriter][HDF5Ext][Material][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Tet::testWriteCellField", "[DataWriter][HDF5Ext][Material][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5ExtMaterial::Hex::testOpenClose", "[DataWriter][HDF5Ext][Material][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Hex::testWriteVertexField", "[DataWriter][HDF5Ext][Material][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5ExtMaterial::Hex::testWriteCellField", "[DataWriter][HDF5Ext][Material][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5ExtMaterial(pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMaterial_Data*
pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMaterial_Data* data = new TestDataWriterHDF5ExtMaterial_Data();assert(data);

    data->opencloseFilename = "tri3_mat.h5";
    data->vertexFilename = "tri3_mat_vertex.h5";
    data->cellFilename = "tri3_mat_cell.h5";

    TestDataWriterMaterial::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMaterial_Data*
pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMaterial_Data* data = new TestDataWriterHDF5ExtMaterial_Data();assert(data);

    data->opencloseFilename = "quad4_mat.h5";
    data->vertexFilename = "quad4_mat_vertex.h5";
    data->cellFilename = "quad4_mat_cell.h5";

    TestDataWriterMaterial::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMaterial_Data*
pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMaterial_Data* data = new TestDataWriterHDF5ExtMaterial_Data();assert(data);

    data->opencloseFilename = "tet4_mat.h5";
    data->vertexFilename = "tet4_mat_vertex.h5";
    data->cellFilename = "tet4_mat_cell.h5";

    TestDataWriterMaterial::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtMaterial_Data*
pylith::meshio::TestDataWriterHDF5ExtMaterial_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtMaterial_Data* data = new TestDataWriterHDF5ExtMaterial_Data();assert(data);

    data->opencloseFilename = "hex8_mat.h5";
    data->vertexFilename = "hex8_mat_vertex.h5";
    data->cellFilename = "hex8_mat_cell.h5";

    TestDataWriterMaterial::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
