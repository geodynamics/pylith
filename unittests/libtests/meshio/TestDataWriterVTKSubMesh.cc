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

#include "TestDataWriterVTKSubMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTKSubMesh );

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;
typedef pylith::topology::Field<pylith::topology::SubMesh> SubMeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKSubMesh::setUp(void)
{ // setUp
  TestDataWriterVTK::setUp();
  _mesh = 0;
  _submesh = 0;
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKSubMesh::tearDown(void)
{ // tearDown
  TestDataWriterVTK::tearDown();
  delete _mesh; _mesh = 0;
  delete _submesh; _submesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTKSubMesh::testConstructor(void)
{ // testConstructor
  DataWriterVTK<topology::SubMesh, MeshField> writer;

  CPPUNIT_ASSERT(0 == writer._viewer);
  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
} // testConstructor

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTKSubMesh::testTimeStep(void)
{ // testTimeStep
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  const double t = _data->time;
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
pylith::meshio::TestDataWriterVTKSubMesh::testWriteVertexField(void)
{ // testWriteVertexField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::SubMesh, MeshField> writer;

  topology::Fields<MeshField> vertexFields(*_mesh);
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numVertexFields;

  const double t = _data->time;
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
pylith::meshio::TestDataWriterVTKSubMesh::testWriteCellField(void)
{ // testWriteCellField
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK<topology::SubMesh, SubMeshField> writer;

  topology::Fields<SubMeshField> cellFields(*_submesh);
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numCellFields;

  const double t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(*_submesh, numTimeSteps);
    writer.openTimeStep(t, *_submesh);
    for (int i=0; i < nfields; ++i) {
      const SubMeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
      writer.writeCellField(t, field);
      CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(*_submesh, numTimeSteps, label, id);
    writer.openTimeStep(t, *_submesh, label, id);
    for (int i=0; i < nfields; ++i) {
      const SubMeshField& field = cellFields.get(_data->cellFieldsInfo[i].name);
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
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKSubMesh::_initialize(void)
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

  CPPUNIT_ASSERT(0 != _data->bcLabel);
  delete _submesh; _submesh = new topology::SubMesh(*_mesh, _data->bcLabel);
  const ALE::Obj<topology::Mesh::SieveMesh>& sieveMesh = _mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<topology::SubMesh::SieveMesh>& sieveSubMesh =
    _submesh->sieveMesh();
  assert(!sieveSubMesh.isNull());
  sieveSubMesh->setRealSection("coordinates", 
			       sieveMesh->getRealSection("coordinates"));
  //_mesh->view("BC mesh");
} // _initialize

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterVTKSubMesh::_createVertexFields(
	    topology::Fields<MeshField>* fields) const
{ // _createVertexFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numVertexFields;

    const ALE::Obj<topology::SubMesh::SieveMesh>& sieveMesh = 
      _mesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    const ALE::Obj<topology::SubMesh::SieveMesh::label_sequence>& vertices =
      sieveMesh->depthStratum(0);
    CPPUNIT_ASSERT(!vertices.isNull());
    const topology::SubMesh::SieveMesh::label_sequence::iterator verticesEnd =
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

      const ALE::Obj<topology::SubMesh::RealSection>& section = field.section();
      CPPUNIT_ASSERT(!section.isNull());
      int ipt = 0;
      for (topology::SubMesh::SieveMesh::label_sequence::iterator v_iter=vertices->begin();
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
pylith::meshio::TestDataWriterVTKSubMesh::_createCellFields(
	     topology::Fields<SubMeshField>* fields) const
{ // _createCellFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numCellFields;

    const ALE::Obj<topology::SubMesh::SieveMesh>& sieveSubMesh =
      _submesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveSubMesh.isNull());
    const ALE::Obj<topology::SubMesh::SieveMesh::label_sequence>& cells = 
      sieveSubMesh->heightStratum(1);
    assert(!cells.isNull());
    const topology::SubMesh::SieveMesh::label_sequence::iterator cellsBegin = 
      cells->begin();
    const topology::SubMesh::SieveMesh::label_sequence::iterator cellsEnd = 
      cells->end();

    // Set cell fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->cellFieldsInfo[i].name;
      const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      SubMeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim, 1);
      field.allocate();
      field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

      const ALE::Obj<topology::SubMesh::RealSection>& section = field.section();
      CPPUNIT_ASSERT(!section.isNull());
      int icell = 0;
      for (topology::SubMesh::SieveMesh::label_sequence::iterator c_iter=cellsBegin;
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
