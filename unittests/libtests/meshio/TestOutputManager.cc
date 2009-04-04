// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestOutputManager.hh" // Implementation of class methods

#include "pylith/meshio/OutputManager.hh"

#include "TestDataWriterVTK.hh" // USES TestDataWriterVTK::checkFile()

#include "pylith/meshio/CellFilterAvg.hh" // USES CellFilterAvg
#include "pylith/meshio/VertexFilterVecNorm.hh" // USES VertexFilterVecNorm
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/feassemble/Quadrature2D.hh" // USES Quadrature

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <math.h> // USES sqrt()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestOutputManager );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestOutputManager::testConstructor(void)
{ // testConstructor
  OutputManager<topology::Mesh> manager;
} // testConstructor

// ----------------------------------------------------------------------
// Test coordsys()
void
pylith::meshio::TestOutputManager::testCoordsys(void)
{ // testCoordsys
  OutputManager<topology::Mesh> manager;

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
  OutputManager<topology::Mesh> manager;

  CPPUNIT_ASSERT(0 == manager._writer);

  DataWriterVTK<topology::Mesh> writer;
  manager.writer(&writer);
  CPPUNIT_ASSERT(0 != manager._writer);
} // testWriter

// ----------------------------------------------------------------------
// Test vertexFilter()
void
pylith::meshio::TestOutputManager::testVertexFilter(void)
{ // testVertexFilter
  OutputManager<topology::Mesh> manager;

  CPPUNIT_ASSERT(0 == manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);

  VertexFilterVecNorm<topology::Mesh> filter;
  manager.vertexFilter(&filter);
  CPPUNIT_ASSERT(0 != manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);
} // testVertexFilter

// ----------------------------------------------------------------------
// Test cellFilter().
void
pylith::meshio::TestOutputManager::testCellFilter(void)
{ // testCellFilter
  OutputManager<topology::Mesh> manager;

  CPPUNIT_ASSERT(0 == manager._vertexFilter);
  CPPUNIT_ASSERT(0 == manager._cellFilter);

  CellFilterAvg<topology::Mesh> filter;
  manager.cellFilter(&filter);
  CPPUNIT_ASSERT(0 != manager._cellFilter);
  CPPUNIT_ASSERT(0 == manager._vertexFilter);
} // testCellFilter

// ----------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestOutputManager::testOpenClose(void)
{ // testOpenClose
  OutputManager<topology::Mesh> manager;

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;

  // TODO Replace DataVTKWriter with writer that has nontrivial
  // open()/close().
  DataWriterVTK<topology::Mesh> writer;
  manager.writer(&writer);

  manager.open(mesh, &cs, numTimeSteps);
  manager.close();
} // testOpenClose

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep().
void
pylith::meshio::TestOutputManager::testOpenCloseTimeStep(void)
{ // testOpenCloseTimeStep
  OutputManager<topology::Mesh> manager;

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const double t = 1.2;
  const char* filenameRoot = "output.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK<topology::Mesh> writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);
  manager.writer(&writer);
  
  manager.open(mesh, &cs, numTimeSteps);
  manager.openTimeStep(t, mesh, &cs);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);
} // testOpenCloseTimeStep

// ----------------------------------------------------------------------
// Test appendVertexField().
void
pylith::meshio::TestOutputManager::testAppendVertexField(void)
{ // testAppendVertexField
  OutputManager<topology::Mesh> manager;

  const char* meshFilename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int nvertices = 4;
  const char* label = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = 
    topology::FieldBase::VECTOR;
  const double fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
    3.1, 3.2,
    4.1, 4.2
  };

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename(meshFilename);
  iohandler.read(&mesh);

  // Set vertex field
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  ALE::Obj<real_section_type> field = 
    new real_section_type(mesh->comm(), mesh->debug());
  field->setChart(mesh->getSieve()->getChart());
  field->setFiberDimension(vertices, fiberDim);
  mesh->allocate(field);

  CPPUNIT_ASSERT_EQUAL(nvertices, int(vertices->size()));

  int ipt = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++ipt) {
    const double* values = &fieldValues[ipt*fiberDim];
    field->updatePoint(*v_iter, values);
  } // for

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const double t = 1.2;
  const char* filenameRoot = "output_vertex.vtk";
  const char* filenameRootF = "output_vertex_filter.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);
  manager.writer(&writer);
  
  manager.open(mesh, &cs, numTimeSteps);
  manager.openTimeStep(t, mesh, &cs);
  manager.appendVertexField(t, fieldName, field, fieldType, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);

  VertexFilterVecNorm filter;
  manager.vertexFilter(&filter);
  writer.filename(filenameRootF);
  manager.writer(&writer);
  
  manager.open(mesh, &cs, numTimeSteps);
  manager.openTimeStep(t, mesh, &cs);
  manager.appendVertexField(t, fieldName, field, fieldType, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);
} // testAppendVertexField

// ----------------------------------------------------------------------
// Test appendCellField().
void
pylith::meshio::TestOutputManager::testAppendCellField(void)
{ // testAppendCellField
  OutputManager manager;

  const char* meshFilename = "data/tri3.mesh";
  const int fiberDim = 2;
  const int ncells = 2;
  const char* fieldName = "field data";
  const topology::FieldBase::VectorFieldEnum fieldType = OTHER_FIELD;
  const double fieldValues[] = {
    1.1, 1.2,
    2.1, 2.2,
  };

  MeshIOAscii iohandler;
  ALE::Obj<Mesh> mesh;
  iohandler.filename(meshFilename);
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

  // Set cell field
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  ALE::Obj<real_section_type> field = 
    new real_section_type(mesh->comm(), mesh->debug());
  field->setChart(mesh->getSieve()->getChart());
  field->setFiberDimension(cells, fiberDim);
  mesh->allocate(field);

  CPPUNIT_ASSERT_EQUAL(ncells, int(cells->size()));

  int ipt = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++ipt) {
    const double* values = &fieldValues[ipt*fiberDim];
    field->updatePoint(*c_iter, values);
  } // for

  spatialdata::geocoords::CSCart cs;
  const int numTimeSteps = 1;
  const double t = 1.2;
  const char* filenameRoot = "output_cell.vtk";
  const char* filenameRootF = "output_cell_filter.vtk";
  const char* timeFormat = "%3.1f";

  DataWriterVTK writer;
  writer.filename(filenameRoot);
  writer.timeFormat(timeFormat);
  manager.writer(&writer);

  manager.open(mesh, &cs, numTimeSteps);
  manager.openTimeStep(t, mesh, &cs);
  manager.appendCellField(t, fieldName, field, fieldType, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRoot, t, timeFormat);

  const int cellDim = 2;
  const int numBasis = 4;
  const int numQuadPts = 2;
  const int spaceDim = 2;
  const double basis[] = {
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
    1.0, 1.0,
  };
  const double basisDerivRef[] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
  };
  const double quadPtsRef[] = {
    1.0, 0.0,
   -1.0, 0.0,};
  const double quadWts[] = { 1.5, 0.5 };
  const double minJacobian = 1.0;

  feassemble::Quadrature2D quadrature;
  quadrature.initialize(basis, basisDerivRef, quadPtsRef, quadWts,
			cellDim, numBasis, numQuadPts, spaceDim);


  CellFilterAvg filter;
  filter.quadrature(&quadrature);
  manager.cellFilter(&filter);
  writer.filename(filenameRootF);
  manager.writer(&writer);
  
  manager.open(mesh, &cs, numTimeSteps);
  manager.openTimeStep(t, mesh, &cs);
  manager.appendCellField(t, fieldName, field, fieldType, mesh);
  manager.closeTimeStep();
  manager.close();

  TestDataWriterVTK::checkFile(filenameRootF, t, timeFormat);
} // testAppendCellField


// End of file 
