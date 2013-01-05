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
// Copyright (c) 2010-2012 University of California, Davis
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

    DM dmMesh = _faultMesh->dmMesh();
    PetscInt       vStart, vEnd;
    PetscErrorCode err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

    // Set vertex fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->vertexFieldsInfo[i].name;
      const int fiberDim = _data->vertexFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      MeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
      field.allocate();
      field.vectorFieldType(_data->vertexFieldsInfo[i].field_type);

      PetscSection section = field.petscSection();
      Vec          vec     = field.localVector();
      PetscScalar *a;
      CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
      err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
      for(PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof, off;

        err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < dof; ++d) {
          a[off+d] = _data->vertexFields[i][(v-vStart)*dof+d];
        }
      } // for
      err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(_data->numVertices, vEnd-vStart);
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

    DM dmMesh = _faultMesh->dmMesh();
    PetscInt       cStart, cEnd;
    const PetscInt height = 0;
    PetscErrorCode err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetHeightStratum(dmMesh, height, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

    // Set cell fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->cellFieldsInfo[i].name;
      const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      MeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
      field.allocate();
      field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

      PetscSection section = field.petscSection();
      Vec          vec     = field.localVector();
      PetscScalar *a;
      CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
      err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
      for(PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt dof, off;

        err = PetscSectionGetDof(section, c, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(section, c, &off);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < dof; ++d) {
          a[off+d] = _data->cellFields[i][(c-cStart)*dof+d];
        }
      } // for
      err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(_data->numCells, cEnd-cStart);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createCellFields


// End of file 
