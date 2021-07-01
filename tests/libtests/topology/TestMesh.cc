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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestMesh);

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestMesh::testConstructor(void) { // testConstructor
    PYLITH_METHOD_BEGIN;

    int result = 0;

    Mesh mesh;
    CPPUNIT_ASSERT(!mesh._dm);
    CPPUNIT_ASSERT_EQUAL(0, mesh.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh.getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);

    int dim = 2;
    Mesh mesh2(dim);
    CPPUNIT_ASSERT(mesh2._dm);
    CPPUNIT_ASSERT_EQUAL(dim, mesh2.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh2.getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    dim = 1;
    Mesh mesh3(dim, PETSC_COMM_SELF);
    CPPUNIT_ASSERT(mesh3._dm);
    CPPUNIT_ASSERT_EQUAL(dim, mesh3.getDimension());
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh3.getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test getDM().
void
pylith::topology::TestMesh::testDMMesh(void) { // testDMMesh
    PYLITH_METHOD_BEGIN;

    const int dim = 2;
    PetscInt dmDim;
    Mesh mesh(dim);

    PetscDM dmMesh = mesh.getDM();CPPUNIT_ASSERT(dmMesh);
    PetscErrorCode err = DMGetDimension(dmMesh, &dmDim);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(dim, dmDim);

    PYLITH_METHOD_END;
} // testDMMesh


// ----------------------------------------------------------------------
// Test getCoordSys().
void
pylith::topology::TestMesh::testCoordsys(void) { // testCoordsys
    PYLITH_METHOD_BEGIN;

    Mesh mesh;

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(2);

    mesh.setCoordSys(&cs);

    CPPUNIT_ASSERT_EQUAL(cs.getSpaceDim(), mesh.getCoordSys()->getSpaceDim());

    PYLITH_METHOD_END;
} // testCoordsys


// ----------------------------------------------------------------------
// Test getDimension().
void
pylith::topology::TestMesh::testDimension(void) { // testDimension
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    CPPUNIT_ASSERT_EQUAL(0, mesh.getDimension());

    const int dim = 2;
    Mesh mesh2(dim);
    CPPUNIT_ASSERT_EQUAL(dim, mesh2.getDimension());

    PYLITH_METHOD_END;
} // testDimension


// ----------------------------------------------------------------------
// Test numCells(), numVertices().
void
pylith::topology::TestMesh::testAccessors(void) { // testAccessors
    PYLITH_METHOD_BEGIN;

    { // Tri
        const char* filename = "data/tri3.mesh";
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.filename(filename);
        iohandler.read(&mesh);

        CPPUNIT_ASSERT_EQUAL(4, pylith::topology::MeshOps::getNumVertices(mesh));
        CPPUNIT_ASSERT_EQUAL(2, pylith::topology::MeshOps::getNumCells(mesh));
    } // Tri

    { // Hex
        const char* filename = "data/twohex8.mesh";
        Mesh mesh;
        meshio::MeshIOAscii iohandler;
        iohandler.filename(filename);
        iohandler.read(&mesh);

        CPPUNIT_ASSERT_EQUAL(12, pylith::topology::MeshOps::getNumVertices(mesh));
        CPPUNIT_ASSERT_EQUAL(2, pylith::topology::MeshOps::getNumCells(mesh));
    } // Hex

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test comm().
void
pylith::topology::TestMesh::testComm(void) { // testComm
    PYLITH_METHOD_BEGIN;

    Mesh mesh;
    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, mesh.getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_IDENT), result);

    Mesh mesh2(2, PETSC_COMM_SELF);
    result = 0;
    MPI_Comm_compare(PETSC_COMM_SELF, mesh2.getComm(), &result);
    CPPUNIT_ASSERT_EQUAL(int(MPI_CONGRUENT), result);

    PYLITH_METHOD_END;
} // testComm


// ----------------------------------------------------------------------
// Test view().
void
pylith::topology::TestMesh::testView(void) { // testView
    PYLITH_METHOD_BEGIN;

    const char* filename = "data/tri3.mesh";

    Mesh mesh;
    meshio::MeshIOAscii iohandler;
    iohandler.filename(filename);
    iohandler.read(&mesh);

    mesh.view();
    mesh.view(":mesh.view:ascii_info_detail");
    mesh.view("vtk:mesh.vtk:ascii_vtk");

    PYLITH_METHOD_END;
} // testView


// End of file
