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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5BCMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterHDF5.hh" // USES DataWriterHDF5
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5BCMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::SubMesh> SubMeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5BCMesh::setUp(void)
{ // setUp
  TestDataWriterBCMesh::setUp();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5BCMesh::tearDown(void)
{ // tearDown
  TestDataWriterBCMesh::tearDown();
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5BCMesh::testConstructor(void)
{ // testConstructor
  DataWriterHDF5<topology::SubMesh, SubMeshField> writer;

  CPPUNIT_ASSERT(0 == writer._viewer);
} // testConstructor

// ----------------------------------------------------------------------
// Test open() and close().
void
pylith::meshio::TestDataWriterHDF5BCMesh::testOpenClose(void)
{ // testOpenClose
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5<topology::SubMesh, SubMeshField> writer;

  writer.filename(_data->timestepFilename);

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_submesh, numTimeSteps);
    writer.openTimeStep(t, *_submesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_submesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_submesh, label, id);
  } // else

  writer.closeTimeStep();
  writer.close();

  checkFile(_data->timestepFilename);
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5BCMesh::testWriteVertexField(void)
{ // testWriteVertexField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5<topology::SubMesh, SubMeshField> writer;

  topology::Fields<SubMeshField> vertexFields(*_submesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);

  const PylithScalr timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numVertexFields;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_submesh, numTimeSteps);
    writer.openTimeStep(t, *_submesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_submesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_submesh, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    SubMeshField& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field, *_submesh);
  } // for
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->vertexFilename);
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5BCMesh::testWriteCellField(void)
{ // testWriteCellField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5<topology::SubMesh, SubMeshField> writer;

  topology::Fields<SubMeshField> cellFields(*_submesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);

  const PylithScalar timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numCellFields;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_submesh, numTimeSteps);
    writer.openTimeStep(t, *_submesh);
    for (int i=0; i < nfields; ++i) {
      SubMeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_submesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_submesh, label, id);
    for (int i=0; i < nfields; ++i) {
      SubMeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field, label, id);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->cellFilename);
} // testWriteCellField


// End of file 
