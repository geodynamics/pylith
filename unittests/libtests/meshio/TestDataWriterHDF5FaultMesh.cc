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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5FaultMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterHDF5.hh" // USES DataWriterHDF5
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <map> // USES std::map

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5FaultMesh );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterFaultMesh::setUp();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5FaultMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  TestDataWriterFaultMesh::tearDown();

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5FaultMesh::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DataWriterHDF5 writer;

  CPPUNIT_ASSERT(!writer._viewer);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterHDF5FaultMesh::testOpenClose(void)
{ // testOpenClose
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterHDF5 writer;

  writer.filename(_data->timestepFilename);

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_faultMesh, numTimeSteps);
    writer.openTimeStep(t, *_faultMesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_faultMesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_faultMesh, label, id);
  } // else

  writer.closeTimeStep();
  writer.close();

  checkFile(_data->timestepFilename);

  PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5FaultMesh::testWriteVertexField(void)
{ // testWriteVertexField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterHDF5 writer;

  topology::Fields vertexFields(*_faultMesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);

  const PylithScalar timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numVertexFields;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_faultMesh, numTimeSteps);
    writer.openTimeStep(t, *_faultMesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_faultMesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_faultMesh, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    topology::Field& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field, *_faultMesh);
  } // for
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->vertexFilename);

  PYLITH_METHOD_END;
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5FaultMesh::testWriteCellField(void)
{ // testWriteCellField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterHDF5 writer;

  topology::Fields cellFields(*_faultMesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);

  const PylithScalar timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numCellFields;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_faultMesh, numTimeSteps);
    writer.openTimeStep(t, *_faultMesh);
    for (int i=0; i < nfields; ++i) {
      topology::Field& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_faultMesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_faultMesh, label, id);
    for (int i=0; i < nfields; ++i) {
      topology::Field& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field, label, id);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->cellFilename);

  PYLITH_METHOD_END;
} // testWriteCellField


// End of file 
