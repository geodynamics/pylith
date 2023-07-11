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

#include "TestDataWriterHDF5Material.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Material_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Material_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5Material_Data* Tri(void);

    static TestDataWriterHDF5Material_Data* Quad(void);

    static TestDataWriterHDF5Material_Data* Tet(void);

    static TestDataWriterHDF5Material_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5Material::Tri::testOpenClose", "[DataWriter][HDF5][Material][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Material::Tri::testWriteVertexField", "[DataWriter][HDF5][Material][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Material::Tri::testWriteCellField", "[DataWriter][HDF5][Material][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Material::Quad::testOpenClose", "[DataWriter][HDF5][Material][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Material::Quad::testWriteVertexField", "[DataWriter][HDF5][Material][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Material::Quad::testWriteCellField", "[DataWriter][HDF5][Material][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Material::Tet::testOpenClose", "[DataWriter][HDF5][Material][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Material::Tet::testWriteVertexField", "[DataWriter][HDF5][Material][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Material::Tet::testWriteCellField", "[DataWriter][HDF5][Material][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterHDF5Material::Hex::testOpenClose", "[DataWriter][HDF5][Material][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5Material::Hex::testWriteVertexField", "[DataWriter][HDF5][Material][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterHDF5Material::Hex::testWriteCellField", "[DataWriter][HDF5][Material][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterHDF5Material(pylith::meshio::TestDataWriterHDF5Material_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Material_Data*
pylith::meshio::TestDataWriterHDF5Material_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Material_Data* data = new TestDataWriterHDF5Material_Data();assert(data);

    data->opencloseFilename = "tri3_mat.h5";
    data->vertexFilename = "tri3_mat_vertex.h5";
    data->cellFilename = "tri3_mat_cell.h5";

    TestDataWriterMaterial::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Material_Data*
pylith::meshio::TestDataWriterHDF5Material_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Material_Data* data = new TestDataWriterHDF5Material_Data();assert(data);

    data->opencloseFilename = "quad4_mat.h5";
    data->vertexFilename = "quad4_mat_vertex.h5";
    data->cellFilename = "quad4_mat_cell.h5";

    TestDataWriterMaterial::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Material_Data*
pylith::meshio::TestDataWriterHDF5Material_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Material_Data* data = new TestDataWriterHDF5Material_Data();assert(data);

    data->opencloseFilename = "tet4_mat.h5";
    data->vertexFilename = "tet4_mat_vertex.h5";
    data->cellFilename = "tet4_mat_cell.h5";

    TestDataWriterMaterial::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5Material_Data*
pylith::meshio::TestDataWriterHDF5Material_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5Material_Data* data = new TestDataWriterHDF5Material_Data();assert(data);

    data->opencloseFilename = "hex8_mat.h5";
    data->vertexFilename = "hex8_mat_vertex.h5";
    data->cellFilename = "hex8_mat_cell.h5";

    TestDataWriterMaterial::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
