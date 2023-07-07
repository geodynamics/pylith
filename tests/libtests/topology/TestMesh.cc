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

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "petscviewerhdf5.h" // USES PetscViewerHDF5

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class TestMesh;
    }
}

class pylith::topology::TestMesh : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Test constructor.
    static
    void testConstructor(void);

    /// Test dmMesh().
    static
    void testDMMesh(void);

    /// Test coordsys().
    static
    void testCoordsys(void);

    /// Test dimension().
    static
    void testDimension(void);

    /// Test numCorners(), numCells(), numVertices(), isSimplex().
    static
    void testAccessors(void);

    /// Test comm().
    static
    void testComm(void);

    /// Test view().
    static
    void testView(void);

}; // class TestMesh

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMesh::testConstructor", "[TestMesh]") {
    pylith::topology::TestMesh::testConstructor();
}
TEST_CASE("TestMesh::testDMMesh", "[TestMesh]") {
    pylith::topology::TestMesh::testDMMesh();
}
TEST_CASE("TestMesh::testCoordsys", "[TestMesh]") {
    pylith::topology::TestMesh::testCoordsys();
}
TEST_CASE("TestMesh::testDimension", "[TestMesh]") {
    pylith::topology::TestMesh::testDimension();
}
TEST_CASE("TestMesh::testAccessors", "[TestMesh]") {
    pylith::topology::TestMesh::testAccessors();
}
TEST_CASE("TestMesh::testComm", "[TestMesh]") {
    pylith::topology::TestMesh::testComm();
}
TEST_CASE("TestMesh::testView", "[TestMesh]") {
    pylith::topology::TestMesh::testView();
}

// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestMesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    int result = 0;

    Mesh mesh;
    CHECK(!mesh._dm);
    CHECK(0 == mesh.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh.getComm(), &result);
    CHECK(int(MPI_IDENT) == result);

    int dim = 2;
    Mesh mesh2(dim);
    CHECK(mesh2._dm);
    CHECK(dim == mesh2.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh2.getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    dim = 1;
    Mesh mesh3(dim, PETSC_COMM_SELF);
    CHECK(mesh3._dm);
    CHECK(dim == mesh3.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh3.getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------

// Test getDM().
void
pylith::topology::TestMesh::testDMMesh(void) { // testDMMesh
    PYLITH_METHOD_BEGIN;

    const int dim = 2;
    PetscInt dmDim;
    Mesh mesh(dim);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscErrorCode err = DMGetDimension(dmMesh, &dmDim);assert(!err);
    CHECK(dim == dmDim);

    PYLITH_METHOD_END;
} // testDMMesh


// ------------------------------------------------------------------------------------------------
// Test getCoordSys().
void
pylith::topology::TestMesh::testCoordsys(void) { // testCoordsys
    PYLITH_METHOD_BEGIN;

    Mesh mesh;

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(2);

    mesh.setCoordSys(&cs);

    CHECK(cs.getSpaceDim() == mesh.getCoordSys()->getSpaceDim());

    PYLITH_METHOD_END;
} // testCoordsys


// ------------------------------------------------------------------------------------------------
// Test getDimension().
void
pylith::topology::TestMesh::testDimension(void) { // testDimension
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    CHECK(0 == mesh.getDimension());

    const int dim = 2;
    Mesh mesh2(dim);
    CHECK(dim == mesh2.getDimension());

    PYLITH_METHOD_END;
} // testDimension


// ------------------------------------------------------------------------------------------------
// Test numCells(), numVertices().
void
pylith::topology::TestMesh::testAccessors(void) { // testAccessors
    PYLITH_METHOD_BEGIN;

    { // Tri
        const char* filename = "data/tri3.mesh";
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.setFilename(filename);
        iohandler.read(&mesh);

        CHECK(4 == pylith::topology::MeshOps::getNumVertices(mesh));
        CHECK(2 == pylith::topology::MeshOps::getNumCells(mesh));
    } // Tri

    { // Hex
        const char* filename = "data/twohex8.mesh";
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.setFilename(filename);
        iohandler.read(&mesh);

        CHECK(12 == pylith::topology::MeshOps::getNumVertices(mesh));
        CHECK(2 == pylith::topology::MeshOps::getNumCells(mesh));
    } // Hex

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestMesh::testComm(void) { // testComm
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh.getComm(), &result);
    CHECK(int(MPI_IDENT) == result);

    Mesh mesh2(2, PETSC_COMM_SELF);
    result = 0;
    MPI_Comm_compare(PETSC_COMM_SELF, mesh2.getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    PYLITH_METHOD_END;
} // testComm


// ------------------------------------------------------------------------------------------------
// Test view().
void
pylith::topology::TestMesh::testView(void) { // testView
    PYLITH_METHOD_BEGIN;

    const char* filename = "data/tri3.mesh";

    Mesh mesh;
    meshio::MeshIOAscii iohandler;
    iohandler.setFilename(filename);
    iohandler.read(&mesh);

    mesh.view();
    mesh.view("ascii:mesh.txt:ascii_info_detail");
    mesh.view("vtk:mesh.vtk:ascii_vtk");
    mesh.view("vtk:mesh.vtu:vtk_vtu");
    mesh.view("ascii:mesh.tex:ascii_latex");
    mesh.view("hdf5:mesh_xdmf.h5:hdf5_xdmf");
    mesh.view("hdf5:mesh_petsc.h5:hdf5_petsc");

    PetscErrorCode err = 0;
    PetscViewer viewer = NULL;
    PetscDM dm = NULL;
    err = PetscViewerHDF5Open(PETSC_COMM_SELF, "mesh_petsc.h5", FILE_MODE_READ, &viewer);REQUIRE(!err);
    err = DMCreate(PETSC_COMM_SELF, &dm);REQUIRE(!err);
    err = DMLoad(dm, viewer);REQUIRE(!err);
    err = DMDestroy(&dm);REQUIRE(!err);
    err = PetscViewerDestroy(&viewer);REQUIRE(!err);

    PYLITH_METHOD_END;
} // testView


// End of file
