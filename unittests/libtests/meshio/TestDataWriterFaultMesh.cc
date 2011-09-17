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

#include "TestDataWriterFaultMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <map> // USES std::map

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::SubMesh> MeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterFaultMesh::setUp(void)
{ // setUp
  _data = 0;
  _mesh = new topology::Mesh();
  _faultMesh = new topology::SubMesh();
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterFaultMesh::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _mesh; _mesh = 0;
  delete _faultMesh; _faultMesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterFaultMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _faultMesh);

  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(_mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->dimension());
  _mesh->coordsys(&cs);

  faults::FaultCohesiveKin fault;
  int firstFaultVertex    = 0;
  int firstLagrangeVertex = _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
  int firstFaultCell      = _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
  const bool constraintCell = true;
  if (constraintCell) {
    firstFaultCell += _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
  }
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(_mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);
  faults::CohesiveTopology::createFaultParallel(_faultMesh, *_mesh, _data->faultId,
						constraintCell);
} // _initialize

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterFaultMesh::_createVertexFields(
	    topology::Fields<MeshField>* fields) const
{ // _createVertexFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _faultMesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numVertexFields;

    const ALE::Obj<topology::SubMesh::SieveMesh>& sieveFaultMesh = 
      _faultMesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveFaultMesh.isNull());
    const ALE::Obj<topology::SubMesh::SieveMesh::label_sequence>& vertices =
      sieveFaultMesh->depthStratum(0);
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
	const PylithScalar* values = &_data->vertexFields[i][ipt*fiberDim];
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
pylith::meshio::TestDataWriterFaultMesh::_createCellFields(
	     topology::Fields<MeshField>* fields) const
{ // _createCellFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numCellFields;

    const ALE::Obj<topology::SubMesh::SieveMesh>& sieveFaultMesh =
      _faultMesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveFaultMesh.isNull());
    const ALE::Obj<topology::SubMesh::SieveMesh::label_sequence>& cells = 
      sieveFaultMesh->heightStratum(0);
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
      MeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
      field.allocate();
      field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

      const ALE::Obj<topology::SubMesh::RealSection>& section = field.section();
      CPPUNIT_ASSERT(!section.isNull());
      int icell = 0;
      for (topology::SubMesh::SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	   c_iter != cellsEnd;
	   ++c_iter, ++icell) {
	const PylithScalar* values = &_data->cellFields[i][icell*fiberDim];
	section->updatePoint(*c_iter, values);
      } // for
      CPPUNIT_ASSERT_EQUAL(_data->numCells, icell);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createCellFields


// End of file 
