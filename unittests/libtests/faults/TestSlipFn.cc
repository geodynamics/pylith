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

#include "TestSlipFn.hh" // Implementation of class methods

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
void
pylith::faults::TestSlipFn::_createFaultMesh(topology::SubMesh* faultMesh,
					     topology::Mesh* mesh,
					     const char* faultLabel,
					     const int faultId)
{ // _createFaultMesh
  PYLITH_METHOD_BEGIN;
  
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(mesh);

  PetscErrorCode err = 0;
  { // Create mesh
    PetscInt firstFaultVertex = 0;
    PetscInt firstLagrangeVertex = 0, firstFaultCell = 0;
    PetscDMLabel groupField = NULL;
    const bool useLagrangeConstraints = true;
    PetscDM faultBoundaryDM = NULL;
    PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    
    err = DMPlexGetStratumSize(dmMesh, faultLabel, 1, &firstLagrangeVertex);CHECK_PETSC_ERROR(err);
    firstFaultCell = firstLagrangeVertex;
    if (useLagrangeConstraints) {
      firstFaultCell += firstLagrangeVertex;
    } // if
    err = DMPlexGetLabel(dmMesh, faultLabel, &groupField);CHECK_PETSC_ERROR(err);
    CPPUNIT_ASSERT(groupField);
    ALE::Obj<SieveFlexMesh> faultBoundary = 0;
    const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
    CPPUNIT_ASSERT(!sieveMesh.isNull());
    CohesiveTopology::createFault(faultMesh, faultBoundary, faultBoundaryDM, *mesh, groupField);
    CohesiveTopology::create(mesh, *faultMesh, faultBoundary, faultBoundaryDM, groupField, faultId, firstFaultVertex, firstLagrangeVertex, firstFaultCell, useLagrangeConstraints);
    err = DMDestroy(&faultBoundaryDM);CHECK_PETSC_ERROR(err);
  } // Create mesh

  { // Need to copy coordinates from mesh to fault mesh since we are not
    // using create() instead of createParallel().
    PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    PetscDM faultDMMesh = faultMesh->dmMesh();CPPUNIT_ASSERT(faultDMMesh);
    const int  spaceDim = mesh->dimension();
    PetscIS subpointIS = NULL;
    const PetscInt *points = NULL;
    PetscSection coordSection = NULL, fcoordSection = NULL;
    PetscInt vStart, vEnd, ffStart, ffEnd;
    
    err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHeightStratum(faultDMMesh, 1, &ffStart, &ffEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);CHECK_PETSC_ERROR(err);
    err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
    err = DMPlexGetCoordinateSection(faultDMMesh, &fcoordSection);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetChart(fcoordSection, vStart, vEnd);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      err = PetscSectionSetDof(fcoordSection, v, spaceDim);CHECK_PETSC_ERROR(err);
    } // for
    err = PetscSectionSetUp(fcoordSection);CHECK_PETSC_ERROR(err);
    PetscVec coordVec, fcoordVec;
    PetscScalar *coords,  *fcoords;
    PetscInt coordSize;
    
    err = PetscSectionGetStorageSize(fcoordSection, &coordSize);CHECK_PETSC_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
    err = VecCreate(mesh->comm(), &fcoordVec);CHECK_PETSC_ERROR(err);
    err = VecSetSizes(fcoordVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
    err = VecSetFromOptions(fcoordVec);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
    err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
    err = VecGetArray(fcoordVec, &fcoords);CHECK_PETSC_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt off, foff;
      
      // Notice that subpointMap[] does not account for cohesive cells
      err = PetscSectionGetOffset(coordSection, points[v]+(ffEnd-ffStart), &off);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(fcoordSection, v, &foff);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < spaceDim; ++d) {
	fcoords[foff+d] = coords[off+d];
      } // for
    } // for
    err = ISRestoreIndices(subpointIS, &points);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(fcoordVec, &fcoords);CHECK_PETSC_ERROR(err);
    err = DMSetCoordinatesLocal(faultDMMesh, fcoordVec);CHECK_PETSC_ERROR(err);
    err = VecDestroy(&fcoordVec);CHECK_PETSC_ERROR(err);
  } // Copy coordiantes


  PYLITH_METHOD_END;
} // _createFaultMesh
