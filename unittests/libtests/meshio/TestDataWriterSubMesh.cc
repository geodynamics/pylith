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

#include "TestDataWriterSubMesh.hh" // Implementation of class methods

#include "data/DataWriterData.hh" // USES DataWriterData

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ----------------------------------------------------------------------
typedef pylith::topology::Field<pylith::topology::Mesh> MeshField;
typedef pylith::topology::Field<pylith::topology::SubMesh> SubMeshField;

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterSubMesh::setUp(void)
{ // setUp
  _data = 0;
  _mesh = 0;
  _submesh = 0;
  _flipFault = false;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterSubMesh::tearDown(void)
{ // tearDown
  delete _data; _data = 0;
  delete _mesh; _mesh = 0;
  delete _submesh; _submesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterSubMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);

  delete _mesh; _mesh = new topology::Mesh;
  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(_mesh);

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->dimension());
  _mesh->coordsys(&cs);

  if (0 != _data->faultLabel) {
    faults::FaultCohesiveKin fault;
    int firstFaultVertex    = 0;
    int firstLagrangeVertex = _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
    int firstFaultCell      = _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
    if (fault.useLagrangeConstraints()) {
      firstFaultCell += _mesh->sieveMesh()->getIntSection(_data->faultLabel)->size();
    }
    fault.label(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(_mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell, _flipFault);
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
pylith::meshio::TestDataWriterSubMesh::_createVertexFields(
	    topology::Fields<MeshField>* fields) const
{ // _createVertexFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numVertexFields;

    DM dmMesh    = _mesh->dmMesh();
    PetscInt       vStart, vEnd;
    PetscErrorCode err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
#if 0
    DM dmSubmesh = _submesh->dmMesh();
    IS subpointIS;
    err = DMPlexCreateSubpointIS(dmSubmesh, &subpointIS);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(subpointIS, &ind);CHECK_PETSC_ERROR(err);
    err = ISRestoreIndices(subpointIS, &ind);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
#endif

    // Set vertex fields
    for (int i=0; i < nfields; ++i) {
      const PetscInt *ind;
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

        err = 
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
  } catch (...) {
    throw;
  } // catch
} // _createVertexFields

// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterSubMesh::_createCellFields(
	     topology::Fields<SubMeshField>* fields) const
{ // _createCellFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(0 != _mesh);
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numCellFields;

    DM              dmMesh = _submesh->dmMesh();
    PetscInt        cStart, cEnd, numCells;
    PetscErrorCode  err;

    CPPUNIT_ASSERT(dmMesh);
    err = DMPlexGetHeightStratum(dmMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    numCells = cEnd - cStart;

    // Set cell fields
    for (int i=0; i < nfields; ++i) {
      const char* name = _data->cellFieldsInfo[i].name;
      const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
      fields->add(name, name);
      SubMeshField& field = fields->get(name);
      field.newSection(topology::FieldBase::CELLS_FIELD, fiberDim, 1);
      field.allocate();
      field.vectorFieldType(_data->cellFieldsInfo[i].field_type);

      PetscSection section = field.petscSection();
      Vec          vec     = field.localVector();
      PetscScalar *a;
      CPPUNIT_ASSERT(section);CPPUNIT_ASSERT(vec);
      err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
      for(PetscInt c = 0; c < numCells; ++c) {
        const PetscInt cell = c+cStart;
        PetscInt       dof, off;

        err = PetscSectionGetDof(section, cell, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(section, cell, &off);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < dof; ++d) {
          a[off+d] = _data->cellFields[i][c*dof+d];
        }
      } // for
      err = VecRestoreArray(vec, &a);CHECK_PETSC_ERROR(err);
      CPPUNIT_ASSERT_EQUAL(_data->numCells, numCells);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } catch (...) {
    throw;
  } // catch
} // _createCellFields


// End of file 
