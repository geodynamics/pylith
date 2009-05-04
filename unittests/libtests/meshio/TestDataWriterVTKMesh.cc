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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKMesh::setUp(void)
{ // setUp
  TestDataWriterVTK::setUp();
  _mesh = 0;
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKMesh::tearDown(void)
{ // tearDown
  TestDataWriterVTK::tearDown();
  delete _mesh; _mesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKMesh::testConstructor(void)
{ // testConstructor
  DataWriterVTK<topology::Mesh, MeshField> writer;

  CPPUNIT_ASSERT(0 == writer._viewer);
  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterVTKMesh::testFilename(void)
{ // testDebug
  DataWriterVTK<topology::Mesh, MeshField> writer;

  const char* filename = "data.vtk";
  writer.filename(filename);
  CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);
} // testFilename

// ----------------------------------------------------------------------
// Test timeFormat()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeFormat(void)
{ // testTimeFormat
  DataWriterVTK<topology::Mesh, MeshField> writer;

  const char* format = "%4.1f";
  writer.timeFormat(format);
  CPPUNIT_ASSERT_EQUAL(std::string(format), writer._timeFormat);
} // testInterpolate

// ----------------------------------------------------------------------
// Test timeConstant()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeConstant(void)
{ // testTimeConstant
  DataWriterVTK<topology::Mesh, MeshField> writer;

  const double value = 4.5;
  writer.timeConstant(value);
  CPPUNIT_ASSERT_EQUAL(value, writer._timeConstant);
} // testInterpolate

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKMesh::testTimeStep(void)
{ // testTimeStep
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::Mesh, MeshField> writer;

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  const double t = _data->time;
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

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  writer.closeTimeStep();
  writer.close();

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  checkFile(_data->timestepFilename, t, _data->timeFormat);
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTKMesh::testWriteVertexField(void)
{ // testWriteVertexField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::Mesh, MeshField> writer;

  topology::Fields<MeshField> vertexFields(*_mesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numVertexFields;

  const double t = _data->time;
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
    const MeshField& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
    writer.writeVertexField(t, field);
    CPPUNIT_ASSERT(writer._wroteVertexHeader);
    CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  } // for
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  
  checkFile(_data->vertexFilename, t, _data->timeFormat);
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTKMesh::testWriteCellField(void)
{ // testWriteCellField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::Mesh, MeshField> writer;

  topology::Fields<MeshField> cellFields(*_mesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numCellFields;

  const double t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_mesh, numTimeSteps);
    writer.openTimeStep(t, *_mesh);
    for (int i=0; i < nfields; ++i) {
      const MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
      CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_mesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_mesh, label, id);
    for (int i=0; i < nfields; ++i) {
      const MeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field, label, id);
      CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  
  checkFile(_data->cellFilename, t, _data->timeFormat);
} // testWriteCellField

// ----------------------------------------------------------------------
// Test _vtkFilename.
void pylith::meshio::TestDataWriterVTKMesh::testVtkFilename(void)
{ // testVtkFilename
  DataWriterVTK<topology::Mesh, MeshField> writer;

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
} // testVtkFilename

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);

  delete _mesh; _mesh = new topology::Mesh;
  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(_mesh);

  if (0 != _data->faultLabel) {
    faults::FaultCohesiveKin fault;
    fault.label(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(_mesh, _flipFault);
  } // if
} // _initialize

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterVTKMesh::_createVertexFields(
	    topology::Fields<MeshField>* fields) const
{ // _createVertexFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numVertexFields;

    const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = _mesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<topology::Mesh::SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const topology::Mesh::SieveMesh::label_sequence::iterator verticesEnd =
      vertices->end();

    // Set vertex fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->vertexFieldsInfo[i].name;
      const int fiberDim = _data->vertexFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      MeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
      field.allocate();
      field.vectorFieldType(_data->vertexFieldsInfo[i].field_type);

      const ALE::Obj<topology::Mesh::RealSection>& section = field.section();
      CPPUNIT_ASSERT(!section.isNull());
      int ipt = 0;
      for (topology::Mesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
	   v_iter != verticesEnd;
	   ++v_iter, ++ipt) {
	const double* values = &_data->vertexFields[i][ipt*fiberDim];
	section->updatePoint(*v_iter, values);
      } // for
      CPPUNIT_ASSERT_EQUAL(_data->numVertices, ipt);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createVertexFields

// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterVTKMesh::_createCellFields(
	     topology::Fields<MeshField>* fields) const
{ // _createCellFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numCellFields;

    const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = _mesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<topology::Mesh::SieveMesh::label_sequence>& cells = 
      (0 == _data->cellsLabel) ? 
      sieveMesh->depthStratum(1) :
      sieveMesh->getLabelStratum(_data->cellsLabel, _data->labelId);
    const topology::Mesh::SieveMesh::label_sequence::iterator cellsEnd = 
      cells->end();

    // Set cell fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->cellFieldsInfo[i].name;
      const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      MeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
      field.allocate();
      field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

      const ALE::Obj<topology::Mesh::RealSection>& section = field.section();
      CPPUNIT_ASSERT(!section.isNull());
      int icell = 0;
      for (topology::Mesh::SieveMesh::label_sequence::iterator c_iter=cells->begin();
	   c_iter != cellsEnd;
	   ++c_iter, ++icell) {
	const double* values = &_data->cellFields[i][icell*fiberDim];
	section->updatePoint(*c_iter, values);
      } // for
      CPPUNIT_ASSERT_EQUAL(_data->numCells, icell);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createCellFields


// End of file 
