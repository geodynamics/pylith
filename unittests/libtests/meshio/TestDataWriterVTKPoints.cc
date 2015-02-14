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

#include "TestDataWriterVTKPoints.hh" // Implementation of class methods

#include "data/DataWriterDataPoints.hh" // USES DataWriterDataPoints

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/OutputSolnPoints.hh" // USES OutputSolnPoints
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKPoints );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKPoints::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  TestDataWriterPoints::setUp();

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKPoints::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  TestDataWriterPoints::tearDown();

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKPoints::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  DataWriterVTK writer;

  CPPUNIT_ASSERT(!writer._viewer);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKPoints::testTimeStep(void)
{ // testTimeStep
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  OutputSolnPoints output;
  DataWriterVTK writer;
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(10.0);

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);
  output.writer(&writer);
  output.setupInterpolator(_mesh, _data->points, _data->numPoints, _data->spaceDim, normalizer);

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    output.open(*_mesh, numTimeSteps);
    output.writePointNames(_data->names, _data->numPoints); // Should to nothing
    output.openTimeStep(t, *_mesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    output.open(*_mesh, numTimeSteps, label, id);
    output.openTimeStep(t, *_mesh, label, id);
  } // else

  output.closeTimeStep();
  output.close();

  // Nothing to check. We do not create VTK files without fields anymore.

  PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKPoints::testWriteVertexField(void)
{ // testWriteVertexField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  OutputSolnPoints output;
  DataWriterVTK writer;
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(10.0);

  topology::Fields vertexFields(*_mesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);
  output.writer(&writer);
  output.setupInterpolator(_mesh, _data->points, _data->numPoints, _data->spaceDim, normalizer);

  const int nfields = _data->numVertexFields;

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (!_data->cellsLabel) {
    output.open(*_mesh, numTimeSteps);
    output.writePointNames(_data->names, _data->numPoints); // Should to nothing
    output.openTimeStep(t, *_mesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    output.open(*_mesh, numTimeSteps, label, id);
    output.openTimeStep(t, *_mesh, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    topology::Field& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    // field.view("FIELD"); // DEBUGGING
    output.appendVertexField(t, field, *_mesh);
    CPPUNIT_ASSERT(writer._wroteVertexHeader);
    CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);
  } // for
  output.closeTimeStep();
  output.close();
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteVertexHeader);
  CPPUNIT_ASSERT_EQUAL(false, writer._wroteCellHeader);

  checkFile(_data->vertexFilename, t, _data->timeFormat);

  PYLITH_METHOD_END;
} // testWriteVertexField


// End of file 
