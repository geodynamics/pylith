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

#include "TestDataWriterVTKMaterial.hh" // Implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMaterial_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKMaterial_Cases {
public:

    // Data factory methods
    static TestDataWriterVTKMaterial_Data* Tri(void);

    static TestDataWriterVTKMaterial_Data* Quad(void);

    static TestDataWriterVTKMaterial_Data* Tet(void);

    static TestDataWriterVTKMaterial_Data* Hex(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestDataWriterVTKMesh::Tri::testTimeStep", "[DataWriter][VTK][Material][Tri][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tri()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMaterial::Tri::testWriteVertexField", "[DataWriter][VTK][Material][Tri][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tri()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMaterial::Tri::testWriteCellField", "[DataWriter][VTK][Material][Tri][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tri()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMaterial::Quad::testTimeStep", "[DataWriter][VTK][Material][Quad][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Quad()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMaterial::Quad::testWriteVertexField", "[DataWriter][VTK][Material][Quad][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Quad()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMaterial::Quad::testWriteCellField", "[DataWriter][VTK][Material][Quad][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Quad()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMaterial::Tet::testTimeStep", "[DataWriter][VTK][Material][Tet][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tet()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMaterial::Tet::testWriteVertexField", "[DataWriter][VTK][Material][Tet][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tet()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMaterial::Tet::testWriteCellField", "[DataWriter][VTK][Material][Tet][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Tet()).testWriteCellField();
}

TEST_CASE("TestDataWriterVTKMaterial::Hex::testTimeStep", "[DataWriter][VTK][Material][Hex][testTimeStep]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Hex()).testTimeStep();
}
TEST_CASE("TestDataWriterVTKMaterial::Hex::testWriteVertexField", "[DataWriter][VTK][Material][Hex][testWriteVertexField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Hex()).testWriteVertexField();
}
TEST_CASE("TestDataWriterVTKMaterial::Hex::testWriteCellField", "[DataWriter][VTK][Material][Hex][testWriteCellField]") {
    pylith::meshio::TestDataWriterVTKMaterial(pylith::meshio::TestDataWriterVTKMaterial_Cases::Hex()).testWriteCellField();
}

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMaterial_Data*
pylith::meshio::TestDataWriterVTKMaterial_Cases::Tri(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMaterial_Data* data = new TestDataWriterVTKMaterial_Data();assert(data);

    data->timestepFilename = "tri3_mat.vtk";
    data->vertexFilename = "tri3_mat_vertex.vtk";
    data->cellFilename = "tri3_mat_cell.vtk";

    TestDataWriterMaterial::setDataTri(data);

    PYLITH_METHOD_RETURN(data);
} // Tri


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMaterial_Data*
pylith::meshio::TestDataWriterVTKMaterial_Cases::Quad(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMaterial_Data* data = new TestDataWriterVTKMaterial_Data();assert(data);

    data->timestepFilename = "quad4_mat.vtk";
    data->vertexFilename = "quad4_mat_vertex.vtk";
    data->cellFilename = "quad4_mat_cell.vtk";

    TestDataWriterMaterial::setDataQuad(data);

    PYLITH_METHOD_RETURN(data);
} // Quad


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMaterial_Data*
pylith::meshio::TestDataWriterVTKMaterial_Cases::Tet(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMaterial_Data* data = new TestDataWriterVTKMaterial_Data();assert(data);

    data->timestepFilename = "tet4_mat.vtk";
    data->vertexFilename = "tet4_mat_vertex.vtk";
    data->cellFilename = "tet4_mat_cell.vtk";

    TestDataWriterMaterial::setDataTet(data);

    PYLITH_METHOD_RETURN(data);
} // Tet


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestDataWriterVTKMaterial_Data*
pylith::meshio::TestDataWriterVTKMaterial_Cases::Hex(void) {
    PYLITH_METHOD_BEGIN;
    TestDataWriterVTKMaterial_Data* data = new TestDataWriterVTKMaterial_Data();assert(data);

    data->timestepFilename = "hex8_mat.vtk";
    data->vertexFilename = "hex8_mat_vertex.vtk";
    data->cellFilename = "hex8_mat_cell.vtk";

    TestDataWriterMaterial::setDataHex(data);

    PYLITH_METHOD_RETURN(data);
} // Hex


// End of file
