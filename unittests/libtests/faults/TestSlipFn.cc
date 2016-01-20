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

#include "TestSlipFn.hh" // Implementation of class methods

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
void
pylith::faults::TestSlipFn::_createFaultMesh(topology::Mesh* faultMesh,
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
    PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    
    err = DMGetStratumSize(dmMesh, faultLabel, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
    firstFaultCell = firstLagrangeVertex;
    if (useLagrangeConstraints) {
      firstFaultCell += firstLagrangeVertex;
    } // if
    err = DMGetLabel(dmMesh, faultLabel, &groupField);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(groupField);
    CohesiveTopology::createFault(faultMesh, *mesh, groupField);
    CohesiveTopology::create(mesh, *faultMesh, groupField, faultId, firstFaultVertex, firstLagrangeVertex, firstFaultCell, useLagrangeConstraints);
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
    
    err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(faultDMMesh, 1, &ffStart, &ffEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinateSection(dmMesh, &coordSection);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinateSection(faultDMMesh, &fcoordSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(fcoordSection, vStart, vEnd);PYLITH_CHECK_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      err = PetscSectionSetDof(fcoordSection, v, spaceDim);PYLITH_CHECK_ERROR(err);
    } // for
    err = PetscSectionSetUp(fcoordSection);PYLITH_CHECK_ERROR(err);
    PetscVec coordVec, fcoordVec;
    PetscScalar *coords,  *fcoords;
    PetscInt coordSize;
    
    err = PetscSectionGetStorageSize(fcoordSection, &coordSize);PYLITH_CHECK_ERROR(err);
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);
    err = VecCreate(mesh->comm(), &fcoordVec);PYLITH_CHECK_ERROR(err);
    err = VecSetSizes(fcoordVec, coordSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
    err = VecSetFromOptions(fcoordVec);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(fcoordVec, &fcoords);PYLITH_CHECK_ERROR(err);
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt off, foff;
      
      // Notice that subpointMap[] does not account for cohesive cells
      err = PetscSectionGetOffset(coordSection, points[v]+(ffEnd-ffStart), &off);PYLITH_CHECK_ERROR(err);
      err = PetscSectionGetOffset(fcoordSection, v, &foff);PYLITH_CHECK_ERROR(err);
      for(PetscInt d = 0; d < spaceDim; ++d) {
	fcoords[foff+d] = coords[off+d];
      } // for
    } // for
    err = ISRestoreIndices(subpointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
    err = VecRestoreArray(coordVec, &coords);PYLITH_CHECK_ERROR(err);
    err = VecRestoreArray(fcoordVec, &fcoords);PYLITH_CHECK_ERROR(err);
    err = DMSetCoordinatesLocal(faultDMMesh, fcoordVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&fcoordVec);PYLITH_CHECK_ERROR(err);
  } // Copy coordiantes


  PYLITH_METHOD_END;
} // _createFaultMesh
