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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMesh );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterMesh::setUp();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  TestDataWriterMesh::tearDown();

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKMesh::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  CPPUNIT_ASSERT(!writer._viewer);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterVTKMesh::testFilename(void)
{ // testDebug
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  const char* filename = "data.vtk";
  writer.filename(filename);
  CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);

  PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test timeFormat()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeFormat(void)
{ // testTimeFormat
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  const char* format = "%4.1f";
  writer.timeFormat(format);
  CPPUNIT_ASSERT_EQUAL(std::string(format), writer._timeFormat);

  PYLITH_METHOD_END;
} // testTimeFormat

// ----------------------------------------------------------------------
// Test timeConstant()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeConstant(void)
{ // testTimeConstant
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  const PylithScalar value = 4.5;
  writer.timeConstant(value);
  CPPUNIT_ASSERT_EQUAL(value, writer._timeConstant);

  PYLITH_METHOD_END;
} // testTimeConstant

// ----------------------------------------------------------------------
// Test precision()
void
pylith::meshio::TestDataWriterVTKMesh::testPrecision(void)
{ // testPrecision
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  const int value = 4;
  writer.precision(value);
  CPPUNIT_ASSERT_EQUAL(value, writer._precision);

  PYLITH_METHOD_END;
} // testPrecision

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeStep(void)
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK writer;

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);

  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
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
pylith::meshio::TestDataWriterVTKMesh::testWriteVertexField(void)
{ // testWriteVertexField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK writer;

  topology::Fields vertexFields(*_mesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numVertexFields;

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    topology::Field& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field, *_mesh);

    // Make sure we can reuse field
    std::string fieldLabel = std::string(field.label()) + std::string("2");
    field.label(fieldLabel.c_str());
    field.dimensionalizeOkay(true);
    field.scale(2.0);
    field.dimensionalize();
    writer.writeVertexField(t, field, *_mesh);

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
pylith::meshio::TestDataWriterVTKMesh::testWriteCellField(void)
{ // testWriteCellField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  DataWriterVTK writer;

  topology::Fields cellFields(*_mesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numCellFields;

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
    for (int i=0; i < nfields; ++i) {
      topology::Field& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
      CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
    for (int i=0; i < nfields; ++i) {
      topology::Field& field = cellFields.get(_data->cellFieldsInfo[i].name);
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

// ----------------------------------------------------------------------
// Test _vtkFilename.
void pylith::meshio::TestDataWriterVTKMesh::testVtkFilename(void)
{ // testVtkFilename
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  // Append info to filename if number of time steps is 0.
  writer._numTimeSteps = 0;
  writer._filename = "output.vtk";
  CPPUNIT_ASSERT_EQUAL(std::string("output_info.vtk"), writer._vtkFilename(0.0));
		       
  // Use default normalization of 1.0, remove period from time stamp.
  writer._numTimeSteps = 100;
  writer._filename = "output.vtk";
  writer.timeFormat("%05.2f");
  CPPUNIT_ASSERT_EQUAL(std::string("output_t0230.vtk"), 
		       writer._vtkFilename(2.3));
  
  // Use normalization of 20.0, remove period from time stamp.
  writer._numTimeSteps = 100;
  writer._filename = "output.vtk";
  writer.timeFormat("%05.2f");
  writer.timeConstant(20.0);
  CPPUNIT_ASSERT_EQUAL(std::string("output_t0250.vtk"), 
		       writer._vtkFilename(50.0));

  PYLITH_METHOD_END;
} // testVtkFilename


// End of file 
