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

#include "TestOutputManager.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputManager.hh"

#include "TestDataWriterVTK.hh" // USES TestDataWriterVTK::checkFile()

#include "pylith/meshio/CellFilterAvg.hh" // USES CellFilterAvg
#include "pylith/meshio/VertexFilterVecNorm.hh" // USES VertexFilterVecNorm
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <math.h> // USES sqrt()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::meshio::TestOutputManager);

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputManager::testConstructor(void) { // testConstructor
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test getCoordSys()
void
pylith::meshio::TestOutputManager::testCoordsys(void) { // testCoordsys
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    CPPUNIT_ASSERT(!manager._coordsys);

    spatialdata::geocoords::CSCart cs;
    manager.setCoordSys(&cs);
    CPPUNIT_ASSERT(0 != manager._coordsys);

    PYLITH_METHOD_END;
} // testCoordsys


// ----------------------------------------------------------------------
// Test writer()
void
pylith::meshio::TestOutputManager::testWriter(void) { // testWriter
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    CPPUNIT_ASSERT(!manager._writer);

    DataWriterVTK writer;
    manager.writer(&writer);
    CPPUNIT_ASSERT(manager._writer);

    PYLITH_METHOD_END;
} // testWriter


// ----------------------------------------------------------------------
// Test vertexFilter()
void
pylith::meshio::TestOutputManager::testVertexFilter(void) { // testVertexFilter
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    CPPUNIT_ASSERT(!manager._vertexFilter);
    CPPUNIT_ASSERT(!manager._cellFilter);

    VertexFilterVecNorm filter;
    manager.vertexFilter(&filter);
    CPPUNIT_ASSERT(manager._vertexFilter);
    CPPUNIT_ASSERT(!manager._cellFilter);

    PYLITH_METHOD_END;
} // testVertexFilter


// ----------------------------------------------------------------------
// Test cellFilter().
void
pylith::meshio::TestOutputManager::testCellFilter(void) { // testCellFilter
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    CPPUNIT_ASSERT(!manager._vertexFilter);
    CPPUNIT_ASSERT(!manager._cellFilter);

    CellFilterAvg filter;
    manager.cellFilter(&filter);
    CPPUNIT_ASSERT(manager._cellFilter);
    CPPUNIT_ASSERT(!manager._vertexFilter);

    PYLITH_METHOD_END;
} // testCellFilter


// ----------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestOutputManager::testOpenClose(void) { // testOpenClose
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    topology::Mesh mesh;
    MeshIOAscii iohandler;
    iohandler.filename("data/tri3.mesh");
    iohandler.read(&mesh);

    spatialdata::geocoords::CSCart cs;
    const bool isInfo = false;

    // TODO Replace DataVTKWriter with writer that has nontrivial
    // open()/close().
    DataWriterVTK writer;
    manager.writer(&writer);

    manager.open(mesh, isInfo);
    manager.close();

    PYLITH_METHOD_END;
} // testOpenClose


// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep().
void
pylith::meshio::TestOutputManager::testOpenCloseTimeStep(void) { // testOpenCloseTimeStep
    PYLITH_METHOD_BEGIN;

    OutputManager manager;

    topology::Mesh mesh;
    MeshIOAscii iohandler;
    iohandler.filename("data/tri3.mesh");
    iohandler.read(&mesh);

    spatialdata::geocoords::CSCart cs;
    const bool isInfo = false;
    const PylithScalar t = 1.2;
    const char* filenameRoot = "output.vtk";
    const char* timeFormat = "%3.1f";

    DataWriterVTK writer;
    writer.filename(filenameRoot);
    writer.timeFormat(timeFormat);
    manager.writer(&writer);

    manager.open(mesh, isInfo);
    manager.openTimeStep(t, mesh);
    manager.closeTimeStep();
    manager.close();

    // Nothing to check. We do not create VTK files without fields anymore.

    PYLITH_METHOD_END;
} // testOpenCloseTimeStep


// ----------------------------------------------------------------------
// Test appendVertexField().
void
pylith::meshio::TestOutputManager::testAppendVertexField(void) { // testAppendVertexField
    PYLITH_METHOD_BEGIN;

    const char* meshFilename = "data/tri3.mesh";
    const int fiberDim = 2;
    const int nvertices = 4;
    const char* label = "field data";
    const topology::FieldBase::VectorFieldEnum fieldType =
        topology::FieldBase::VECTOR;
    const PylithScalar fieldValues[] = {
        1.1, 1.2,
        2.1, 2.2,
        3.1, 3.2,
        4.1, 4.2
    };
    const PylithScalar scale = 2.0;

    topology::Mesh mesh;
    MeshIOAscii iohandler;
    iohandler.filename(meshFilename);
    iohandler.read(&mesh);

    // Set vertex field
    PetscDM dmMesh = mesh.getDM();CPPUNIT_ASSERT(dmMesh);
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    topology::Field field(mesh);
    field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    field.allocate();
    field.setLabel(label);
    field.vectorFieldType(fieldType);
    field.scale(scale);

    topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

    for (PetscInt v = vStart, index = 0; v < vEnd; ++v) {
        const PetscInt off = fieldVisitor.sectionOffset(v);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
        for (PetscInt d = 0; d < fiberDim; ++d, ++index) {
            fieldArray[off+d] = fieldValues[index]/scale;
        } // for
    } // for
    CPPUNIT_ASSERT_EQUAL(nvertices, vEnd-vStart);

    spatialdata::geocoords::CSCart cs;
    const bool isInfo = false;
    const PylithScalar t = 1.2;
    const char* filenameRoot = "output_vertex.vtk";
    const char* filenameRootF = "output_vertex_filter.vtk";
    const char* timeFormat = "%3.1f";

    DataWriterVTK writer;
    writer.filename(filenameRoot);
    writer.timeFormat(timeFormat);

    OutputManager manager;
    manager.writer(&writer);
    manager.open(mesh, isInfo);
    manager.openTimeStep(t, mesh);
    manager.appendVertexField(t, field, mesh);
    manager.closeTimeStep();
    manager.close();

    TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);

    VertexFilterVecNorm filter;
    manager.vertexFilter(&filter);
    writer.filename(filenameRootF);
    manager.writer(&writer);

    manager.open(mesh, isInfo);
    manager.openTimeStep(t, mesh);
    manager.appendVertexField(t, field, mesh);
    manager.closeTimeStep();
    manager.close();

    TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);

    PYLITH_METHOD_END;
} // testAppendVertexField


// ----------------------------------------------------------------------
// Test appendCellField().
void
pylith::meshio::TestOutputManager::testAppendCellField(void) { // testAppendCellField
    PYLITH_METHOD_BEGIN;

    const char* meshFilename = "data/tri3.mesh";
    const int fiberDim = 2;
    const int ncells = 2;
    const char* label = "field data";
    const topology::FieldBase::VectorFieldEnum fieldType =
        topology::FieldBase::MULTI_SCALAR;
    const PylithScalar fieldValues[] = {
        1.1, 1.2,
        2.1, 2.2,
    };
    const PylithScalar scale = 4.0;

    topology::Mesh mesh;
    MeshIOAscii iohandler;
    iohandler.filename(meshFilename);
    iohandler.read(&mesh);

    // Set cell field
    PetscDM dmMesh = mesh.getDM();CPPUNIT_ASSERT(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    PetscInt numCells = cellsStratum.size();

    topology::Field field(mesh);
    field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
    field.allocate();
    field.setLabel(label);
    field.vectorFieldType(fieldType);
    field.scale(scale);

    topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

    for (PetscInt c = 0, index = 0; c < numCells; ++c) {
        const PetscInt cell = c+cStart;

        const PetscInt off = fieldVisitor.sectionOffset(cell);
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(cell));
        for (PetscInt d = 0; d < fiberDim; ++d, ++index) {
            fieldArray[off+d] = fieldValues[index]/scale;
        } // for
    } // for
    CPPUNIT_ASSERT_EQUAL(ncells, cEnd-cStart);

    spatialdata::geocoords::CSCart cs;
    const bool isInfo = false;
    const PylithScalar t = 1.2;
    const char* filenameRoot = "output_cell.vtk";
    const char* filenameRootF = "output_cell_filter.vtk";
    const char* timeFormat = "%3.1f";

    DataWriterVTK writer;
    writer.filename(filenameRoot);
    writer.timeFormat(timeFormat);

    OutputManager manager;
    manager.writer(&writer);
    manager.open(mesh, isInfo);
    manager.openTimeStep(t, mesh);
    manager.appendCellField(t, field);
    manager.closeTimeStep();
    manager.close();

    TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);

    const int cellDim = 2;
    const int numBasis = 4;
    const int numQuadPts = 2;
    const int spaceDim = 2;
    const PylithScalar basis[] = {
        1.0, 1.0,
        1.0, 1.0,
        1.0, 1.0,
        1.0, 1.0,
    };
    const PylithScalar basisDerivRef[] = {
        1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0,
    };
    const PylithScalar quadPtsRef[] = {
        1.0, 0.0,
        -1.0, 0.0,
    };
    const PylithScalar quadWts[] = { 1.5, 0.5 };

    feassemble::Quadrature quadrature;
    quadrature.initialize(basis, numQuadPts, numBasis,
                          basisDerivRef, numQuadPts, numBasis, cellDim,
                          quadPtsRef, numQuadPts, cellDim,
                          quadWts, numQuadPts,
                          spaceDim);

    CellFilterAvg filter;
    filter.quadrature(&quadrature);
    manager.cellFilter(&filter);
    writer.filename(filenameRootF);
    manager.writer(&writer);

    manager.open(mesh, isInfo);
    manager.openTimeStep(t, mesh);
    manager.appendCellField(t, field);
    manager.closeTimeStep();
    manager.close();

    TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);

    PYLITH_METHOD_END;
} // testAppendCellField


// End of file
