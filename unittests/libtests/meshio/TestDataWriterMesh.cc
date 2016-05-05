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

#include "TestDataWriterMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterMesh::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _mesh = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterMesh::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _mesh; _mesh = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterMesh::_initialize(void)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  delete _mesh; _mesh = new topology::Mesh;CPPUNIT_ASSERT(_mesh);
  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(_mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->dimension());
  _mesh->coordsys(&cs);

  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(10.0);
  topology::MeshOps::nondimensionalize(_mesh, normalizer);

  if (_data->faultLabel) {
    faults::FaultCohesiveKin fault;
    const bool useLagrangeConstraints = true;
    PetscInt firstFaultVertex = 0;
    PetscInt firstLagrangeVertex = 0, firstFaultCell = 0;
    PetscErrorCode err = DMGetStratumSize(_mesh->dmMesh(), _data->faultLabel, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
    firstFaultCell = firstLagrangeVertex;
    if (useLagrangeConstraints) {
      firstFaultCell += firstLagrangeVertex;
    } // if
    fault.label(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(_mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  } // if

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterMesh::_createVertexFields(topology::Fields* fields) const
{ // _createVertexFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  const int nfields = _data->numVertexFields;
  
  PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);  
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();
  
  // Set vertex fields
  for (int i=0; i < nfields; ++i) {
    const char* name = _data->vertexFieldsInfo[i].name;
    const int fiberDim = _data->vertexFieldsInfo[i].fiber_dim;
    fields->add(name, name);
    topology::Field& field = fields->get(name);
    field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    field.allocate();
    field.vectorFieldType(_data->vertexFieldsInfo[i].field_type);

    topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
    
    for(PetscInt v = vStart, index=0; v < vEnd; ++v) {
      const PetscInt off = fieldVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
      for(PetscInt d = 0; d < fiberDim; ++d, ++index) {
	fieldArray[off+d] = _data->vertexFields[i][index];
      } // for
    } // for
    CPPUNIT_ASSERT_EQUAL(_data->numVertices, vEnd-vStart);
  } // for

  PYLITH_METHOD_END;
} // _createVertexFields

// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterMesh::_createCellFields(topology::Fields* fields) const
{ // _createCellFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);

  const int nfields = _data->numCellFields;

  PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);  
  topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();
  PetscInt numCells = cellsStratum.size();

  topology::StratumIS* cellsIS = (_data->cellsLabel) ? new topology::StratumIS(dmMesh, _data->cellsLabel, _data->labelId) : 0;
  const PetscInt *cells  = NULL;
  if (cellsIS) {
    numCells = cellsIS->size();
    cells = cellsIS->points();
  } // if

  // Set cell fields
  for (int i=0; i < nfields; ++i) {
    const char* name = _data->cellFieldsInfo[i].name;
    const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
    fields->add(name, name);
    topology::Field& field = fields->get(name);
    field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
    field.allocate();
    field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

    topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);
    
    for(PetscInt c = 0, index = 0; c < numCells; ++c) {
      const PetscInt cell = cells ? cells[c] : c+cStart;
      
      const PetscInt off = fieldVisitor.sectionOffset(cell);
      CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(cell));
      for(PetscInt d = 0; d < fiberDim; ++d, ++index) {
	fieldArray[off+d] = _data->cellFields[i][index];
      } // for
    } // for
    CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);
  } // for
  delete cellsIS; cellsIS = 0;

  PYLITH_METHOD_END;
} // _createCellFields


// End of file 
