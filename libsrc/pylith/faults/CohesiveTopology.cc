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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "CohesiveTopology.hh" // implementation of object methods

#include "TopologyOps.hh" // USES TopologyOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include <cassert> // USES assert()

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

extern "C" PetscErrorCode DMPlexMarkBoundaryFaces_Internal(DM, PetscInt, DMLabel);

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::createFault(topology::Mesh* faultMesh,
                                              const topology::Mesh& mesh,
                                              PetscDMLabel groupField)
{ // createFault
  PYLITH_METHOD_BEGIN;

  assert(faultMesh);
  PetscErrorCode err;

  faultMesh->coordsys(mesh.coordsys());
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);

  PetscInt dim, depth, gdepth;
  err = DMPlexGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);

  // Convert fault to a DM
  err = MPI_Allreduce(&depth, &gdepth, 1, MPIU_INT, MPI_MAX, mesh.comm());PYLITH_CHECK_ERROR(err);
  if (gdepth == dim) {
    PetscDM subdm = NULL;
    PetscDMLabel label = NULL;
    const char *groupName = "", *labelName = "boundary";

    if (groupField) {err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);}
    err = DMPlexCreateSubmesh(dmMesh, groupField, 1, &subdm);PYLITH_CHECK_ERROR(err);
    err = DMPlexOrient(subdm);PYLITH_CHECK_ERROR(err);

    err = DMPlexCreateLabel(subdm, labelName);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabel(subdm, labelName, &label);PYLITH_CHECK_ERROR(err);
    err = DMPlexMarkBoundaryFaces_Internal(subdm, 1, label);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(subdm, label);PYLITH_CHECK_ERROR(err);
    std::string submeshLabel = "fault_" + std::string(groupName);
    faultMesh->dmMesh(subdm, submeshLabel.c_str());
  } else {
    PetscDM faultDMMeshTmp = NULL, faultDMMesh = NULL;
    PetscDMLabel subpointMapTmp = NULL, subpointMap = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    PetscInt depth, newDepth, h, numPoints = 0, p;
    const char *groupName = "";

    if (groupField) {err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);}
    err = DMPlexCreateSubmesh(dmMesh, groupField, 1, &faultDMMeshTmp);PYLITH_CHECK_ERROR(err);
    err = DMPlexInterpolate(faultDMMeshTmp, &faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetVTKCellHeight(faultDMMeshTmp, &h);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetVTKCellHeight(faultDMMesh, h);PYLITH_CHECK_ERROR(err);
    err = DMPlexOrient(faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexCopyCoordinates(faultDMMeshTmp, faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(faultDMMeshTmp, &subpointMapTmp);PYLITH_CHECK_ERROR(err);
    err = DMLabelCreate("subpoint_map", &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(subpointMapTmp, 0, &pointIS);PYLITH_CHECK_ERROR(err);
    if (pointIS) {
      err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
      err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    }
    for (p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(subpointMap, points[p], 0);PYLITH_CHECK_ERROR(err);
    }
    if (pointIS) {err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);}
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepth(faultDMMeshTmp, &depth);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepth(faultDMMesh, &newDepth);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(subpointMapTmp, depth, &pointIS);PYLITH_CHECK_ERROR(err);
    if (pointIS) {
      err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
      err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    }
    for (p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(subpointMap, points[p], newDepth);PYLITH_CHECK_ERROR(err);
    }
    if (pointIS) {err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);}
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetSubpointMap(faultDMMesh, subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelDestroy(&subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&faultDMMeshTmp);PYLITH_CHECK_ERROR(err);
    
    std::string submeshLabel = "fault_" + std::string(groupName);
    faultMesh->dmMesh(faultDMMesh, submeshLabel.c_str());

    PetscDMLabel label = NULL;
    const char *labelName = "boundary";

    err = DMPlexCreateLabel(faultDMMesh, labelName);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabel(faultDMMesh, labelName, &label);PYLITH_CHECK_ERROR(err);
    err = DMPlexMarkBoundaryFaces(faultDMMesh, label);PYLITH_CHECK_ERROR(err);
  }

  PYLITH_METHOD_END;
} // createFault

void
pylith::faults::CohesiveTopology::create(topology::Mesh* mesh,
						     const topology::Mesh& faultMesh,
                             PetscDMLabel faultBdLabel,
						     const int materialId,
						     int& firstFaultVertex,
						     int& firstLagrangeVertex,
						     int& firstFaultCell,
						     const bool constraintCell)
{ // createInterpolated
  assert(mesh);
  PetscDM        sdm = NULL;
  PetscDM        dm  = mesh->dmMesh();assert(dm);
  PetscDMLabel   subpointMap = NULL, label = NULL, mlabel = NULL;
  PetscInt       dim, cMax, cEnd, numCohesiveCellsOld;
  PetscErrorCode err;

  // Have to remember the old number of cohesive cells
  err = DMPlexGetHeightStratum(dm, 0, NULL, &cEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetHybridBounds(dm, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
  numCohesiveCellsOld = cEnd - (cMax < 0 ? cEnd : cMax);
  // Create cohesive cells
  err = DMPlexGetSubpointMap(faultMesh.dmMesh(), &subpointMap);PYLITH_CHECK_ERROR(err);
  err = DMLabelDuplicate(subpointMap, &label);PYLITH_CHECK_ERROR(err);
  err = DMLabelClearStratum(label, mesh->dimension());PYLITH_CHECK_ERROR(err);
  // Completes the set of cells scheduled to be replaced
  err = DMPlexLabelCohesiveComplete(dm, label, faultBdLabel, PETSC_FALSE, faultMesh.dmMesh());PYLITH_CHECK_ERROR(err);
  err = DMPlexConstructCohesiveCells(dm, label, &sdm);PYLITH_CHECK_ERROR(err);

  err = DMPlexGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetLabel(sdm, "material-id", &mlabel);PYLITH_CHECK_ERROR(err);
  if (mlabel) {
    err = DMPlexGetHeightStratum(sdm, 0, NULL, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(sdm, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
    assert(cEnd > cMax + numCohesiveCellsOld);
    for (PetscInt cell = cMax; cell < cEnd - numCohesiveCellsOld; ++cell) {
      PetscInt onBd;

      /* Eliminate hybrid cells on the boundary of the split from cohesive label,
         they are marked with -(cell number) since the hybrid cell number aliases vertices in the old mesh */
      err = DMLabelGetValue(label, -cell, &onBd);PYLITH_CHECK_ERROR(err);
      //if (onBd == dim) continue;
      err = DMLabelSetValue(mlabel, cell, materialId);PYLITH_CHECK_ERROR(err);
    }
  }
  err = DMLabelDestroy(&label);PYLITH_CHECK_ERROR(err);

  PetscReal lengthScale = 1.0;
  err = DMPlexGetScale(dm, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
  err = DMPlexSetScale(sdm, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(sdm);
} // createInterpolated

// ----------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::CohesiveTopology::createFaultParallel(topology::Mesh* faultMesh,
						      const topology::Mesh& mesh,
						      const int materialId,
						      const char* label,
						      const bool constraintCell)
{ // createFaultParallel
  PYLITH_METHOD_BEGIN;

  assert(faultMesh);
  const char    *labelname = "material-id";
  PetscErrorCode err;

  faultMesh->coordsys(mesh.coordsys());

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscDM dmFaultMesh;


  err = DMPlexCreateCohesiveSubmesh(dmMesh, constraintCell ? PETSC_TRUE : PETSC_FALSE, labelname, materialId, &dmFaultMesh);PYLITH_CHECK_ERROR(err);
  err = DMPlexOrient(dmFaultMesh);PYLITH_CHECK_ERROR(err);
  std::string meshLabel = "fault_" + std::string(label);

  PetscReal lengthScale = 1.0;
  err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
  err = DMPlexSetScale(dmFaultMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

  faultMesh->dmMesh(dmFaultMesh, meshLabel.c_str());

  PYLITH_METHOD_END;
} // createFaultParallel


// End of file
