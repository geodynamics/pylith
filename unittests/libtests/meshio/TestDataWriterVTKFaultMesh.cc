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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterVTKFaultMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <map> // USES std::map

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKFaultMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::SubMesh> MeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterFaultMesh::setUp();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKFaultMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  TestDataWriterFaultMesh::tearDown();

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKFaultMesh::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  CPPUNIT_ASSERT(!writer._viewer);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKFaultMesh::testTimeStep(void)
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);

  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

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

  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  writer.closeTimeStep();
  writer.close();

  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  // Nothing to check. We do not create VTK files without fields anymore.

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKFaultMesh::testWriteVertexField(void)
{ // testWriteVertexField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  topology::Fields<MeshField> vertexFields(*_faultMesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numVertexFields;

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
  for (int i=0; i < nfields; ++i) {
    MeshField& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field, *_faultMesh);
    CPPUNIT_ASSERT(writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
  } // for
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
  
  checkFile(_data->vertexFilename, t, _data->timeFormat);

  PYLITH_METHOD_END;
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTKFaultMesh::testWriteCellField(void)
{ // testWriteCellField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  topology::Fields<MeshField> cellFields(*_faultMesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numCellFields;

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_faultMesh, numTimeSteps);
    writer.openTimeStep(t, *_faultMesh);
    for (int i=0; i < nfields; ++i) {
      MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
      CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_faultMesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_faultMesh, label, id);
    for (int i=0; i < nfields; ++i) {
      MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field, label, id);
      CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
  
  checkFile(_data->cellFilename, t, _data->timeFormat);

  PYLITH_METHOD_END;
} // testWriteCellField


// End of file 
