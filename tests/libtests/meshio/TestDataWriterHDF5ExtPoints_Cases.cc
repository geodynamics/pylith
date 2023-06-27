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

#include "TestDataWriterHDF5ExtPoints.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtPoints_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtPoints_Cases {
public:

    // Data factory methods
    static TestDataWriterHDF5ExtPoints_Data* Tri(void);

    static TestDataWriterHDF5ExtPoints_Data* Quad(void);

    static TestDataWriterHDF5ExtPoints_Data* Tet(void);

    static TestDataWriterHDF5ExtPoints_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterHDF5ExtPoints::Tri::testOpenClose", "[DataWriter][HDF5Ext][Points][Tri][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tri()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtPoints::Tri::testWriteVertexField", "[DataWriter][HDF5Ext][Points][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tri()).testWriteVertexField();
}

TEST_CASE("TestDataWriterHDF5ExtPoints::Quad::testOpenClose", "[DataWriter][HDF5Ext][Points][Quad][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Quad()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtPoints::Quad::testWriteVertexField", "[DataWriter][HDF5Ext][Points][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Quad()).testWriteVertexField();
}

TEST_CASE("TestDataWriterHDF5ExtPoints::Tet::testOpenClose", "[DataWriter][HDF5Ext][Points][Tet][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tet()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtPoints::Tet::testWriteVertexField", "[DataWriter][HDF5Ext][Points][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tet()).testWriteVertexField();
}

TEST_CASE("TestDataWriterHDF5ExtPoints::Hex::testOpenClose", "[DataWriter][HDF5Ext][Points][Hex][testOpenClose]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Hex()).testOpenClose();
}
TEST_CASE("TestDataWriterHDF5ExtPoints::Hex::testWriteVertexField", "[DataWriter][HDF5Ext][Points][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterHDF5ExtPoints(pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Hex()).testWriteVertexField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtPoints_Data*
pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtPoints_Data* data = new TestDataWriterHDF5ExtPoints_Data();assert(data);

    data->opencloseFilename = "tri3_points.h5";
    data->vertexFilename = "tri3_points_vertex.h5";

    TestDataWriterPoints::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtPoints_Data*
pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtPoints_Data* data = new TestDataWriterHDF5ExtPoints_Data();assert(data);

    data->opencloseFilename = "quad4_points.h5";
    data->vertexFilename = "quad4_points_vertex.h5";

    TestDataWriterPoints::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtPoints_Data*
pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtPoints_Data* data = new TestDataWriterHDF5ExtPoints_Data();assert(data);

    data->opencloseFilename = "tet4_points.h5";
    data->vertexFilename = "tet4_points_vertex.h5";

    TestDataWriterPoints::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterHDF5ExtPoints_Data*
pylith::meshio::TestDataWriterHDF5ExtPoints_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterHDF5ExtPoints_Data* data = new TestDataWriterHDF5ExtPoints_Data();assert(data);

    data->opencloseFilename = "hex8_points.h5";
    data->vertexFilename = "hex8_points_vertex.h5";

    TestDataWriterPoints::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
