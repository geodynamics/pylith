// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestOutputManager.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
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
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputManager );

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputManager::testConstructor(void)
{ // testConstructor
  OutputManager<topology::Mesh, MeshField> manager;
} // testConstructor

// ----------------------------------------------------------------------
// Test coordsys()
void
pylith::meshio::TestOutputManager::testCoordsys(void)
{ // testCoordsys
  OutputManager<topology::Mesh, MeshField> manager;

  CPPUNIT_ASSERT(0 == manager._coordsys);

  spatialdata::geocoords::CSCart cs;
  manager.coordsys(&cs);
  CPPUNIT_ASSERT(0 != manager._coordsys);
} // testCoordsys

// ----------------------------------------------------------------------
// Test writer()
void
pylith::meshio::TestOutputManager::testWriter(void)
{ // testWriter
  OutputManager<topology::Mesh, MeshField> manager;

  CPPUNIT_ASSERT(0 == manager._writer);

  DataWriterVTK<topology::Mesh, MeshField> writer;
  manager.writer(&writer);
  CPPUNIT_ASSERT(0 != manager._writer);
} // testWriter

// ----------------------------------------------------------------------
// Test vertexFilter()
void
pylith::meshio::TestOutputManager::testVertexFilter(void)
{ // testVertexFilter
  OutputManager<topology::Mesh, MeshField> manager;

  CPPUNIT_ASSERT(0 == manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);

  VertexFilterVecNorm<MeshField> filter;
  manager.vertexFilter(&filter);
  CPPUNIT_ASSERT(0 != manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);
} // testVertexFilter

// ----------------------------------------------------------------------
// Test cellFilter().
void
pylith::meshio::TestOutputManager::testCellFilter(void)
{ // testCellFilter
  OutputManager<topology::Mesh, MeshField> manager;

  CPPUNIT_ASSERT(0 == manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);

  CellFilterAvg<topology::Mesh, MeshField> filter;
  manager.cellFilter(&filter);
  CPPUNIT_ASSERT(0 != manager._cellFilter);
  CPPUNIT_ASSERT(0 == manager._vertexFilter);
} // testCellFilter

// ----------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestOutputManager::testOpenClose(void)
{ // testOpenClose
  OutputManager<topology::Mesh, MeshField> manager;

  topology::Mesh mesh;
  MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;

  // TODO Replace DataVTKWriter with writer that has nontrivial
  // open()/close().
  DataWriterVTK<topology::Mesh, MeshField> writer;
  manager.writer(&writer);

  manager.open(mesh, numTimeSteps);
  manager.close();
} // testOpenClose

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep().
void
pylith::meshio::TestOutputManager::testOpenCloseTimeStep(void)
{ // testOpenCloseTimeStep
  OutputManager<topology::Mesh, MeshField> manager;

  topology::Mesh mesh;
  MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const PylithScalar t = 1.2;
  const char* filenameRoot = "output.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK<topology::Mesh, MeshField> writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);
  manager.writer(&writer);
  
  manager.open(mesh, numTimeSteps);
  manager.openTimeStep(t, mesh);
  manager.closeTimeStep();
  manager.close();

  // We do not create empty VTK files anymore
  //TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);
} // testOpenCloseTimeStep

// ----------------------------------------------------------------------
// Test appendVertexField().
void
pylith::meshio::TestOutputManager::testAppendVertexField(void)
{ // testAppendVertexField
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
  DM dmMesh = mesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  MeshField field(mesh);
  field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  field.allocate();
  field.label(label);
  field.vectorFieldType(fieldType);
  field.scale(scale);
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  PetscScalar *a;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off;

    err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      a[off+d] = fieldValues[(v-vStart)*dof+d]/scale;
    }
  } // for
  err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(nvertices, vEnd-vStart);

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const PylithScalar t = 1.2;
  const char* filenameRoot = "output_vertex.vtk";
  const char* filenameRootF = "output_vertex_filter.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK<topology::Mesh, MeshField> writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);

  OutputManager<topology::Mesh, MeshField> manager;
  manager.writer(&writer);
  manager.open(mesh, numTimeSteps);
  manager.openTimeStep(t, mesh);
  manager.appendVertexField(t, field, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);

  VertexFilterVecNorm<MeshField> filter;
  manager.vertexFilter(&filter);
  writer.filename(filenameRootF);
  manager.writer(&writer);
  
  manager.open(mesh, numTimeSteps);
  manager.openTimeStep(t, mesh);
  manager.appendVertexField(t, field, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);
} // testAppendVertexField

// ----------------------------------------------------------------------
// Test appendCellField().
void
pylith::meshio::TestOutputManager::testAppendCellField(void)
{ // testAppendCellField

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
  DM              dmMesh = mesh.dmMesh();
  PetscInt        cStart, cEnd, numCells;
  PetscErrorCode  err;

  CPPUNIT_ASSERT(dmMesh);
  err = DMComplexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  numCells = cEnd - cStart;

  MeshField field(mesh);
  field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  field.allocate();
  field.label(label);
  field.vectorFieldType(fieldType);
  field.scale(scale);
  PetscSection section = field.petscSection();
  Vec          vec     = field.localVector();
  PetscScalar *a;
  CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
  err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = c+cStart;
    PetscInt       dof, off;

    err = PetscSectionGetDof(section, cell, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, cell, &off);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      a[off+d] = fieldValues[c*dof+d]/scale;
    }
  } // for
  err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(ncells, cEnd-cStart);

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const PylithScalar t = 1.2;
  const char* filenameRoot = "output_cell.vtk";
  const char* filenameRootF = "output_cell_filter.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK<topology::Mesh, MeshField> writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);

  OutputManager<topology::Mesh, MeshField> manager;
  manager.writer(&writer);
  manager.open(mesh, numTimeSteps);
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
   -1.0, 0.0,};
  const PylithScalar quadWts[] = { 1.5, 0.5 };
  const PylithScalar minJacobian = 1.0;

  feassemble::Quadrature<topology::Mesh> quadrature;
  quadrature.initialize(basis, numQuadPts, numBasis, 
			basisDerivRef, numQuadPts, numBasis, cellDim,
			quadPtsRef, numQuadPts, cellDim,
			quadWts, numQuadPts,
			spaceDim);

  CellFilterAvg<topology::Mesh, MeshField> filter;
  filter.quadrature(&quadrature);
  manager.cellFilter(&filter);
  writer.filename(filenameRootF);
  manager.writer(&writer);
  
  manager.open(mesh, numTimeSteps);
  manager.openTimeStep(t, mesh);
  manager.appendCellField(t, field);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);
} // testAppendCellField


// End of file 
