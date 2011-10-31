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

#include "TestDataWriterHDF5ExtMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include "pylith/utils/types.hh" // HASA PylithScalar

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5ExtMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::setUp(void)
{ // setUp
  TestDataWriterMesh::setUp();
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::tearDown(void)
{ // tearDown
  TestDataWriterMesh::tearDown();
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testConstructor(void)
{ // testConstructor
  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  CPPUNIT_ASSERT(writer._h5);
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testFilename(void)
{ // testDebug
  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  const char* filename = "data.h5";
  writer.filename(filename);
  CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);
} // testFilename

// ----------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testOpenClose(void)
{ // testOpenClose
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  writer.filename(_data->timestepFilename);

  const PylithScalar t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
  } // else

  writer.close();

  checkFile(_data->timestepFilename);
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteVertexField(void)
{ // testWriteVertexField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  topology::Fields<MeshField> vertexFields(*_mesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);

  const PylithScalar timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numVertexFields;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    MeshField& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field, *_mesh);
  } // for
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->vertexFilename);
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteCellField(void)
{ // testWriteCellField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  topology::Fields<MeshField> cellFields(*_mesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);

  const PylithScalar timeScale = 4.0;
  writer.timeScale(timeScale);
  const PylithScalar t = _data->time / timeScale;

  const int nfields = _data->numCellFields;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
    for (int i=0; i < nfields; ++i) {
      MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
    for (int i=0; i < nfields; ++i) {
      MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field, label, id);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  
  checkFile(_data->cellFilename);
} // testWriteCellField

// ----------------------------------------------------------------------
// Test _hdf5Filename().
void pylith::meshio::TestDataWriterHDF5ExtMesh::testHdf5Filename(void)
{ // testHdf5Filename
  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  // Append info to filename if number of time steps is 0.
  writer._numTimeSteps = 0;
  writer._filename = "output.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_info.h5"), writer._hdf5Filename());
		       
  writer._numTimeSteps = 5;
  writer._filename = "output_abc.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abc.h5"),
		       writer._hdf5Filename());
  
  writer._numTimeSteps = 10;
  writer._filename = "output_abcd.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abcd.h5"), 
		       writer._hdf5Filename());
} // testHdf5Filename

// ----------------------------------------------------------------------
// Test _datasetFilename().
void pylith::meshio::TestDataWriterHDF5ExtMesh::testDatasetFilename(void)
{ // testDatasetFilename
  DataWriterHDF5Ext<topology::Mesh, MeshField> writer;

  // Append info to filename if number of time steps is 0.
  writer._numTimeSteps = 0;
  writer._filename = "output.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_info_ABCD.dat"),
		       writer._datasetFilename("ABCD"));
		       
  writer._numTimeSteps = 5;
  writer._filename = "output_abc.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abc_field1.dat"),
		       writer._datasetFilename("field1"));
  
  writer._numTimeSteps = 10;
  writer._filename = "output_abcd.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abcd_field2.dat"), 
		       writer._datasetFilename("field2"));
} // testDatasetFilename


// End of file 
