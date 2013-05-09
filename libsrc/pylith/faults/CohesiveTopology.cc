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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "CohesiveTopology.hh" // implementation of object methods

#include "TopologyOps.hh" // USES TopologyOps
#include "TopologyVisitors.hh" // USES TopologyVisitors
#include "pylith/topology/SubMesh.hh" // USES SubMesh

#include <Selection.hh> // Algorithms for submeshes

#include <cassert> // USES assert()

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

extern PetscErrorCode DMPlexGetOrientedFace(DM dm, PetscInt cell, PetscInt faceSize, const PetscInt face[], PetscInt numCorners, PetscInt indices[], PetscInt origVertices[], PetscInt faceVertices[], PetscBool *posOriented);

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;
typedef pylith::topology::Mesh::IntSection IntSection;

// Alleviate the need to call the very slooooow mesh->stratify()
// routine by setting the depth/height of new vertices and cells
// manually.
#define FAST_STRATIFY 1

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::createFault(topology::SubMesh* faultMesh,
                                              ALE::Obj<SieveFlexMesh>& faultBoundary,
                                              DM& faultBoundaryDM,
                                              const topology::Mesh& mesh,
                                              DMLabel groupField,
                                              const bool flipFault)
{ // createFault
  PYLITH_METHOD_BEGIN;

  assert(0 != faultMesh);
  assert(groupField);

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("SerialFaultCreation");
  PetscErrorCode err;

  faultMesh->coordsys(mesh);
  DM       dmMesh = mesh.dmMesh();
  PetscInt dim, depth;

  assert(dmMesh);
  err = DMPlexGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);

  // Convert fault to a DM
  if (depth == dim) {
    DM                 subdm;
    DMLabel            label;
    const char        *groupName, *labelName = "boundary";
    std::ostringstream tmp;

    err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubmesh(dmMesh, groupName, 1, &subdm);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateLabel(subdm, labelName);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabel(subdm, labelName, &label);PYLITH_CHECK_ERROR(err);
    err = DMPlexMarkBoundaryFaces(subdm, label);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubmesh(subdm, labelName, 1, &faultBoundaryDM);PYLITH_CHECK_ERROR(err);
    faultMesh->setDMMesh(subdm);
  } else {
    DM              faultDMMeshTmp, faultDMMesh;
    DMLabel         subpointMapTmp, subpointMap;
    IS              pointIS;
    const PetscInt *points;
    PetscInt        depth, newDepth, h, numPoints, p;
    const char     *groupName;

    err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubmesh(dmMesh, groupName, 1, &faultDMMeshTmp);PYLITH_CHECK_ERROR(err);
    err = DMPlexInterpolate(faultDMMeshTmp, &faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetVTKCellHeight(faultDMMeshTmp, &h);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetVTKCellHeight(faultDMMesh, h);PYLITH_CHECK_ERROR(err);
    err = DMPlexOrient(faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexCopyCoordinates(faultDMMeshTmp, faultDMMesh);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(faultDMMeshTmp, &subpointMapTmp);PYLITH_CHECK_ERROR(err);
    err = DMLabelCreate("subpoint_map", &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(subpointMapTmp, 0, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    for (p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(subpointMap, points[p], 0);PYLITH_CHECK_ERROR(err);
    }
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepth(faultDMMeshTmp, &depth);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepth(faultDMMesh, &newDepth);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(subpointMapTmp, depth, &pointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    for (p = 0; p < numPoints; ++p) {
      err = DMLabelSetValue(subpointMap, points[p], newDepth);PYLITH_CHECK_ERROR(err);
    }
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetSubpointMap(faultDMMesh, subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelDestroy(&subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&faultDMMeshTmp);PYLITH_CHECK_ERROR(err);
    if (flipFault) {
      PetscInt maxConeSize, *revcone, *revconeO;
      PetscInt pStart, pEnd;

      err = DMPlexGetHeightStratum(faultDMMesh, h, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
      err = DMPlexGetMaxSizes(faultDMMesh, &maxConeSize, NULL);PYLITH_CHECK_ERROR(err);
      err = DMGetWorkArray(faultDMMesh, maxConeSize, PETSC_INT, &revcone);PYLITH_CHECK_ERROR(err);
      err = DMGetWorkArray(faultDMMesh, maxConeSize, PETSC_INT, &revconeO);PYLITH_CHECK_ERROR(err);
      for (PetscInt p = pStart; p < pEnd; ++p) {
        const PetscInt *cone, *coneO;
        PetscInt        coneSize, faceSize, c;

        err = DMPlexGetConeSize(faultDMMesh, p, &coneSize);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetCone(faultDMMesh, p, &cone);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetConeOrientation(faultDMMesh, p, &coneO);PYLITH_CHECK_ERROR(err);
        for (c = 0; c < coneSize; ++c) {
          err = DMPlexGetConeSize(faultDMMesh, cone[coneSize-1-c], &faceSize);PYLITH_CHECK_ERROR(err);
          revcone[c]  = cone[coneSize-1-c];
          revconeO[c] = coneO[coneSize-1-c] >= 0 ? -(faceSize-coneO[coneSize-1-c]) : faceSize+coneO[coneSize-1-c];
        }
        err = DMPlexSetCone(faultDMMesh, p, revcone);PYLITH_CHECK_ERROR(err);
        err = DMPlexSetConeOrientation(faultDMMesh, p, revconeO);PYLITH_CHECK_ERROR(err);
      }
      err = DMRestoreWorkArray(faultDMMesh, maxConeSize, PETSC_INT, &revcone);PYLITH_CHECK_ERROR(err);
      err = DMRestoreWorkArray(faultDMMesh, maxConeSize, PETSC_INT, &revconeO);PYLITH_CHECK_ERROR(err);
    }
    faultMesh->setDMMesh(faultDMMesh);

    DMLabel            label;
    const char        *labelName = "boundary";

    err = DMPlexCreateLabel(faultDMMesh, labelName);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabel(faultDMMesh, labelName, &label);PYLITH_CHECK_ERROR(err);
    err = DMPlexMarkBoundaryFaces(faultDMMesh, label);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubmesh(faultDMMesh, labelName, 1, &faultBoundaryDM);PYLITH_CHECK_ERROR(err);
  }

  logger.stagePop();

  PYLITH_METHOD_END;
} // createFault

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::create(topology::Mesh* mesh,
                                         const topology::SubMesh& faultMesh,
                                         const ALE::Obj<SieveFlexMesh>& faultBoundary,
                                         DM faultBoundaryDM,
                                         DMLabel groupField,
                                         const int materialId,
                                         int& firstFaultVertex,
                                         int& firstLagrangeVertex,
                                         int& firstFaultCell,
                                         const bool constraintCell)
{ // create
  PYLITH_METHOD_BEGIN;

  assert(0 != mesh);
  assert(groupField);

  typedef ALE::SieveAlg<SieveFlexMesh> sieveAlg;
  typedef ALE::Selection<SieveFlexMesh> selection;
  const char    *groupName;
  PetscMPIInt    rank;
  PetscErrorCode err;

  err = MPI_Comm_rank(mesh->comm(), &rank);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);
  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("SerialFaultCreation");

  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  /* DMPlex */
  DM complexMesh = mesh->dmMesh();
  DM faultDMMesh = faultMesh.dmMesh();
  assert(complexMesh);assert(faultDMMesh);

  const int  depth = sieveMesh->depth();
  assert(!sieveMesh->heightStratum(0).isNull());
  const int numCells = sieveMesh->heightStratum(0)->size();
  PetscInt cStart, cEnd;
#define USE_DMCOMPLEX_ON
#ifdef USE_DMCOMPLEX_ON
  err = DMPlexGetHeightStratum(complexMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
  assert(numCells == cEnd - cStart);
#endif
  int numCorners = 0; // The number of vertices in a mesh cell
  int faceSize = 0; // The number of vertices in a mesh face
  PetscInt faceSizeDM;
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int* indices = 0; // The indices of a face vertex set in a cell
  PetscInt *indicesDM = PETSC_NULL; // The indices of a face vertex set in a cell
  const int debug = mesh->debug();
  int oppositeVertex = 0;    // For simplices, the vertex opposite a given face
  TopologyOps::PointArray origVertices;
  TopologyOps::PointArray faceVertices;
  TopologyOps::PointArray neighborVertices;
  PetscInt *origVerticesDM;
  PetscInt *faceVerticesDM;

  if (!rank) {
    if (!sieveMesh.isNull()) {
      numCorners = sieveMesh->getNumCellCorners();
      faceSize   = selection::numFaceVertices(sieveMesh);
      indices    = new int[faceSize];
    }
    PetscInt cellDim, numCornersDM;

    err = DMPlexGetDimension(complexMesh, &cellDim);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(complexMesh, cStart, &numCornersDM);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetNumFaceVertices(complexMesh, cellDim, numCornersDM, &faceSizeDM);PYLITH_CHECK_ERROR(err);
    if (!sieveMesh.isNull()) {
      assert(numCorners == numCornersDM);
      assert(faceSize == faceSizeDM);
    }
    err = PetscMalloc(faceSizeDM * sizeof(PetscInt), &indicesDM);PYLITH_CHECK_ERROR(err);
    /* TODO: Do we need faceSize at all? Blaise was using a nice criterion */
    PetscInt fStart, numFaultCornersDM;

    err = DMPlexGetHeightStratum(faultDMMesh, 1, &fStart, NULL);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(faultDMMesh, fStart, &numFaultCornersDM);PYLITH_CHECK_ERROR(err);
    numFaultCorners = numFaultCornersDM;
  }

  // Add new shadow vertices and possibly Lagrange multipler vertices
#ifdef USE_DMCOMPLEX_ON
  IS              fVertexIS;
  const PetscInt *fVerticesDM;
  PetscInt        numFaultVerticesDM, vStart, vEnd;

  err = DMLabelGetStratumIS(groupField, 1, &fVertexIS);PYLITH_CHECK_ERROR(err);
  err = ISGetLocalSize(fVertexIS, &numFaultVerticesDM);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(fVertexIS, &fVerticesDM);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetDepthStratum(complexMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
#endif
  const ALE::Obj<std::set<std::string> >& groupNames = sieveMesh->getIntSections();
  assert(!groupNames.isNull());
  std::map<point_type,point_type> vertexRenumber;
  std::map<point_type,point_type> vertexLagrangeRenumber;
  std::map<point_type,point_type> cellRenumber;
  std::map<point_type,point_type> vertexRenumberDM;
  std::map<point_type,point_type> vertexLagrangeRenumberDM;
  std::map<point_type,point_type> cellRenumberDM;
  PetscInt ffStart, ffEnd, numFaultFacesDM, fvtStart, fvtEnd;
  PetscInt faultVertexOffsetDM, firstFaultVertexDM, firstLagrangeVertexDM, firstFaultCellDM, extraVertices, extraCells;

  err = DMPlexGetDepthStratum(faultDMMesh, 0, &fvtStart, &fvtEnd);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetHeightStratum(faultDMMesh, 1, &ffStart, &ffEnd);PYLITH_CHECK_ERROR(err);
  numFaultFacesDM = ffEnd - ffStart;
#ifdef USE_DMCOMPLEX_ON
  /* TODO This will have to change for multiple faults */
  PetscInt numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices;

  extraVertices         = numFaultVerticesDM * (constraintCell ? 2 : 1); /* Total number of fault vertices on this fault (shadow + Lagrange) */
  extraCells            = numFaultFacesDM;                               /* Total number of fault cells */
  firstFaultVertexDM    = vEnd + extraCells;
  firstLagrangeVertexDM = firstFaultVertexDM + firstLagrangeVertex;
  firstFaultCellDM      = cEnd;
  mesh->getPointTypeSizes(&numNormalCells, &numCohesiveCells, &numNormalVertices, &numShadowVertices, &numLagrangeVertices);
  faultVertexOffsetDM   = numCohesiveCells;
  if (!numNormalCells) {
    mesh->setPointTypeSizes(cEnd-cStart, extraCells, vEnd-vStart, firstLagrangeVertex, constraintCell ? firstLagrangeVertex : 0);
  } else {
    mesh->setPointTypeSizes(numNormalCells, numCohesiveCells+extraCells, numNormalVertices, numShadowVertices+firstLagrangeVertex, constraintCell ? numLagrangeVertices+firstLagrangeVertex : 0);
  }
#endif
  if (firstFaultVertex == 0) {
    firstFaultVertex    += sieve->getBaseSize() + sieve->getCapSize();
    firstLagrangeVertex += firstFaultVertex;
    firstFaultCell      += firstFaultVertex;

#ifdef USE_DMCOMPLEX_ON
    PetscInt pStart, pEnd;
    err = DMPlexGetChart(complexMesh, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    assert(firstFaultVertex == pEnd - pStart);
#endif
  }

#ifdef USE_DMCOMPLEX_ON
  /* DMPlex */
  DM        newMesh;
  PetscInt *newCone;
  PetscInt  maxConeSize = 0;

  err = DMCreate(sieveMesh->comm(), &newMesh);PYLITH_CHECK_ERROR(err);
  err = DMSetType(newMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
  err = DMPlexSetDimension(newMesh, sieveMesh->getDimension());PYLITH_CHECK_ERROR(err);
  err = DMPlexSetChart(newMesh, 0, firstFaultVertexDM + extraVertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt coneSize;
    err = DMPlexGetConeSize(complexMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetConeSize(newMesh, c, coneSize);PYLITH_CHECK_ERROR(err);
    maxConeSize = PetscMax(maxConeSize, coneSize);
  }
  for(PetscInt c = cEnd; c < cEnd+numFaultFacesDM; ++c) {
    err = DMPlexSetConeSize(newMesh, c, constraintCell ? faceSizeDM*3 : faceSizeDM*2);PYLITH_CHECK_ERROR(err);
  }
  err = DMSetUp(newMesh);PYLITH_CHECK_ERROR(err);
  err = PetscMalloc(maxConeSize * sizeof(PetscInt), &newCone);PYLITH_CHECK_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone;
    PetscInt        coneSize, cp;

    err = DMPlexGetCone(complexMesh, c, &cone);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(complexMesh, c, &coneSize);PYLITH_CHECK_ERROR(err);
    for(cp = 0; cp < coneSize; ++cp) {
      newCone[cp] = cone[cp] + extraCells;
    }
    err = DMPlexSetCone(newMesh, c, newCone);PYLITH_CHECK_ERROR(err);
  }
  err = PetscFree(newCone);PYLITH_CHECK_ERROR(err);
  PetscInt cMax, vMax;

  err = DMPlexGetHybridBounds(newMesh, &cMax, PETSC_NULL, PETSC_NULL, &vMax);PYLITH_CHECK_ERROR(err);
  if (cMax < 0) {
    err = DMPlexSetHybridBounds(newMesh, firstFaultCellDM, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
  }
  err = DMPlexSetHybridBounds(newMesh, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DETERMINE, firstLagrangeVertexDM);PYLITH_CHECK_ERROR(err);

  // Renumber labels
  std::set<std::string> names(groupNames->begin(), groupNames->end());
  names.insert(names.begin(), "material-id");
  for(std::set<std::string>::const_iterator name = names.begin(); name != names.end(); ++name) {
    const char     *lname = (*name).c_str();
    IS              idIS;
    PetscInt        n;
    PetscBool       hasLabel;
    const PetscInt *ids;

    err = DMPlexHasLabel(complexMesh, lname, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) continue;
    err = DMPlexGetLabelSize(complexMesh, lname, &n);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetLabelIdIS(complexMesh, lname, &idIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(idIS, &ids);PYLITH_CHECK_ERROR(err);
    for(PetscInt i = 0; i < n; ++i) {
      const PetscInt  id = ids[i];
      const PetscInt *points;
      IS              sIS;
      PetscInt        size;

      err = DMPlexGetStratumSize(complexMesh, lname, id, &size);PYLITH_CHECK_ERROR(err);
      err = DMPlexGetStratumIS(complexMesh, lname, id, &sIS);PYLITH_CHECK_ERROR(err);
      err = ISGetIndices(sIS, &points);PYLITH_CHECK_ERROR(err);
      for(PetscInt s = 0; s < size; ++s) {
        if ((points[s] >= vStart) && (points[s] < vEnd)) {
          err = DMPlexSetLabelValue(newMesh, lname, points[s]+extraCells, id);PYLITH_CHECK_ERROR(err);
        } else {
          err = DMPlexSetLabelValue(newMesh, lname, points[s], id);PYLITH_CHECK_ERROR(err);
        }
      }
      err = ISRestoreIndices(sIS, &points);PYLITH_CHECK_ERROR(err);
      err = ISDestroy(&sIS);PYLITH_CHECK_ERROR(err);
    }
    err = ISRestoreIndices(idIS, &ids);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&idIS);PYLITH_CHECK_ERROR(err);
  } // for
#endif

  // Add fault vertices to groups and construct renumberings
  PetscInt    numLabels;
  std::string skipA = "depth";
  std::string skipB = "material-id";

  err = DMPlexGetNumLabels(complexMesh, &numLabels);PYLITH_CHECK_ERROR(err);
  for (PetscInt fv = 0; fv < numFaultVerticesDM; ++fv, ++firstFaultVertexDM) {
    const PetscInt v    = fVerticesDM[fv];
    const PetscInt vnew = v+extraCells;

    vertexRenumberDM[vnew] = firstFaultVertexDM;
    err = DMPlexSetLabelValue(newMesh, groupName, firstFaultVertexDM, 1);PYLITH_CHECK_ERROR(err);
    if (constraintCell) {
      vertexLagrangeRenumberDM[vnew] = firstLagrangeVertexDM;
      err = DMPlexSetLabelValue(newMesh, groupName, firstLagrangeVertexDM, 1);PYLITH_CHECK_ERROR(err);
      ++firstLagrangeVertexDM;
    } // if

    // Add shadow vertices to other groups, don't add constraint
    // vertices (if they exist) because we don't want BC, etc to act
    // on constraint vertices
    for (PetscInt l = 0; l < numLabels; ++l) {
      const char *name;
      PetscInt    value;

      err = DMPlexGetLabelName(complexMesh, l, &name);PYLITH_CHECK_ERROR(err);
      if (std::string(name) == skipA) continue;
      if (std::string(name) == skipB) continue;

      err = DMPlexGetLabelValue(newMesh, name, vnew, &value);PYLITH_CHECK_ERROR(err);
      if (value != -1) {
        err = DMPlexSetLabelValue(newMesh, name, vertexRenumberDM[vnew], value);PYLITH_CHECK_ERROR(err);
      }
    }
  } // for
  for (PetscInt fv = 0; fv < numFaultVerticesDM; ++fv) {
    const PetscInt fvertex = fVerticesDM[fv];

    vertexRenumber[fvertex] = firstFaultVertex;
    if (debug) std::cout << "Duplicating " << fvertex << " to "	<< vertexRenumber[fvertex] << std::endl;

    logger.stagePop();
    logger.stagePush("SerialFaultStratification");
    // Add shadow and constraint vertices (if they exist) to group
    // associated with fault
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultVertex, 1);
    sieveMesh->setDepth(firstFaultVertex, 0);
#endif
    if (constraintCell) {
      vertexLagrangeRenumber[fvertex] = firstLagrangeVertex;
#if defined(FAST_STRATIFY)
      // OPTIMIZATION
      sieveMesh->setHeight(firstLagrangeVertex, 1);
      sieveMesh->setDepth(firstLagrangeVertex, 0);
#endif
      ++firstLagrangeVertex;
    } // if
    logger.stagePop();
    logger.stagePush("SerialFaultCreation");

    // Add shadow vertices to other groups, don't add constraint
    // vertices (if they exist) because we don't want BC, etc to act
    // on constraint vertices
    const std::set<std::string>::const_iterator namesEnd = groupNames->end();
    for(std::set<std::string>::const_iterator name = groupNames->begin(); name != namesEnd;	++name) {
      const ALE::Obj<IntSection>& group = sieveMesh->getIntSection(*name);
      assert(!group.isNull());
      if (group->getFiberDimension(fvertex))
        group->addPoint(firstFaultVertex, 1);
    } // for
  } // for
  logger.stagePop();
  logger.stagePush("SerialFaultCreation");
  const std::set<std::string>::const_iterator namesEnd = groupNames->end();
  for(std::set<std::string>::const_iterator name = groupNames->begin(); name != namesEnd; ++name) {
    sieveMesh->reallocate(sieveMesh->getIntSection(*name));
  } // for
  logger.stagePop();
  logger.stagePush("SerialFault");

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<SieveFlexMesh::label_type>& material = sieveMesh->getLabel("material-id");
  assert(!material.isNull());
  const int firstCohesiveCell = firstFaultCell;
  const PetscInt firstCohesiveCellDM = firstFaultCellDM;
  TopologyOps::PointSet replaceCells;
  TopologyOps::PointSet noReplaceCells;
  TopologyOps::PointSet replaceVertices;
  TopologyOps::PointSet replaceVerticesDM;
  std::set<SieveFlexMesh::point_type> faceSet;
  PetscInt *cohesiveCone;

  ALE::Obj<SieveSubMesh::sieve_type> ifaultSieve;
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> *sV2 = NULL;
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> *cV2 = NULL;
  SieveSubMesh::label_sequence::iterator *f_iter = NULL;
  IS subpointIS;
  const PetscInt *subpointMap;

  err = DMPlexCreateSubpointIS(faultDMMesh, &subpointIS);PYLITH_CHECK_ERROR(err);
  err = ISGetIndices(subpointIS, &subpointMap);PYLITH_CHECK_ERROR(err);
  err = PetscMalloc3(faceSizeDM,PetscInt,&origVerticesDM,faceSizeDM,PetscInt,&faceVerticesDM,faceSizeDM*3,PetscInt,&cohesiveCone);PYLITH_CHECK_ERROR(err);
  for (PetscInt faceDM = ffStart; faceDM < ffEnd; ++faceDM, ++firstFaultCell, ++firstFaultCellDM) {
    PetscInt face = faceDM;
    if (debug) std::cout << "Considering fault face " << faceDM << std::endl;
    if (!ifaultSieve.isNull()) {
      face = **f_iter; ++(*f_iter);
      ifaultSieve->support(face, *sV2);
      const point_type *cells = sV2->getPoints();
      point_type cell = cells[0];
      point_type otherCell = cells[1];
    }
    const PetscInt *support;
    err = DMPlexGetSupport(faultDMMesh, faceDM, &support);PYLITH_CHECK_ERROR(err);
    // Transform to original mesh cells
    point_type cell      = subpointMap[support[0]];
    point_type otherCell = subpointMap[support[1]];

    if (debug) std::cout << "  Checking orientation against cell " << cell << std::endl;
    PetscInt *faceConeDM = NULL, closureSize, coneSizeDM = 0;
    err = DMPlexGetTransitiveClosure(faultDMMesh, faceDM, PETSC_TRUE, &closureSize, &faceConeDM);PYLITH_CHECK_ERROR(err);
    for (PetscInt c = 0; c < closureSize*2; c += 2) {
      if ((faceConeDM[c] >= fvtStart) && (faceConeDM[c] < fvtEnd)) {
        // TODO After everything works, I will reform the subpointMap with new original vertices
        faceConeDM[coneSizeDM++] = subpointMap[faceConeDM[c]];
      }
    }
    int coneSize;
    const point_type *faceCone;
    if (!ifaultSieve.isNull()) {
      ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*ifaultSieve, face, *cV2);
      coneSize = cV2->getSize();
      faceCone = cV2->getPoints();
      assert(coneSize == coneSizeDM);
    } else {
      faceCone = NULL;
    }
    //ifaultSieve->cone(face, cV2);
    //const int coneSize = cV2.getSize() ? cV2.getSize() : 1;
    //const point_type *faceCone = cV2.getSize() ? cV2.getPoints() : &face;
    bool found = true;

    err = DMPlexGetOrientedFace(complexMesh, cell, coneSizeDM, faceConeDM, numCorners, indicesDM, origVerticesDM, faceVerticesDM, PETSC_NULL);PYLITH_CHECK_ERROR(err);
    if (faceCone) {
      for(int i = 0; i < coneSize; ++i) {
        faceSet.insert(faceCone[i]);
        assert(faceConeDM[i] == faceCone[i]+faultVertexOffsetDM);
      }
      selection::getOrientedFace(sieveMesh, cell, &faceSet, numCorners, indices, &origVertices, &faceVertices);
      if (faceVertices.size() != coneSize) {
        std::cout << "Invalid size for faceVertices " << faceVertices.size() << " of face " << face << " should be " << coneSize << std::endl;
        std::cout << "  firstCohesiveCell " << firstCohesiveCell << " firstFaultCell " << firstFaultCell << " numFaces " << (ffEnd-ffStart) << std::endl;
        std::cout << "  faceSet:" << std::endl;
        for(std::set<SieveFlexMesh::point_type>::const_iterator p_iter = faceSet.begin(); p_iter != faceSet.end(); ++p_iter) {
          std::cout << "    " << *p_iter << std::endl;
        } // if
        std::cout << "  cell cone:" << std::endl;
        ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> cV(std::max(1, sieve->getMaxConeSize()));
        sieve->cone(cell, cV);
        const int coneSize2 = cV.getSize();
        const point_type *cellCone  = cV.getPoints();

        for(int c = 0; c < coneSize2; ++c) std::cout << "    " << cellCone[c] << std::endl;
        std::cout << "  fault cell support:" << std::endl;
        ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV(std::max(1, ifaultSieve->getMaxSupportSize()));
        ifaultSieve->support(face, sV);
        const int supportSize2 = sV.getSize();
        const point_type *cellSupport  = sV.getPoints();
        for(int s = 0; s < supportSize2; ++s) std::cout << "    " << cellSupport[s] << std::endl;
      } // if
      assert(faceVertices.size() == coneSize);
      faceSet.clear();
      for(PetscInt c = 0; c < coneSize; ++c) {
        assert(faceVertices[c]+faultVertexOffsetDM == faceVerticesDM[c]);
      }
    }

    if (numFaultCorners == 0) {
      found = false;
    } else if (numFaultCorners == 2) {
      if (faceVerticesDM[0] != faceConeDM[0])
        found = false;
    } else {
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVerticesDM[v] != faceConeDM[0]))
        ++v;
      for(int c = 0; c < coneSizeDM; ++c, ++v) {
        if (debug) std::cout << "    Checking " << faceConeDM[c] << " against " << faceVerticesDM[v%numFaultCorners] << std::endl;
        if (faceVerticesDM[v%numFaultCorners] != faceConeDM[c]) {
          found = false;
          break;
        } // if
      } // for
    } // if/else

    if (found) {
      if (debug) std::cout << "  Choosing other cell" << std::endl;
      point_type tmpCell = otherCell;
      otherCell = cell;
      cell = tmpCell;
    } else {
      if (debug) std::cout << "  Verifing reverse orientation" << std::endl;
      found = true;
      int v = 0;
      if (numFaultCorners > 0) {
        // Locate first vertex
        while((v < numFaultCorners) && (faceVerticesDM[v] != faceConeDM[coneSizeDM-1]))
          ++v;
        for(int c = coneSizeDM-1; c >= 0; --c, ++v) {
          if (debug) std::cout << "    Checking " << faceConeDM[c] << " against " << faceVerticesDM[v%numFaultCorners] << std::endl;
          if (faceVerticesDM[v%numFaultCorners] != faceConeDM[c]) {
            found = false;
            break;
          } // if
        } // for
      } // if
      if (!found) {
        std::cout << "Considering fault face " << faceDM << std::endl;
        std::cout << "  bordered by cells " << cell << " and " << otherCell << std::endl;
        for(int c = 0; c < coneSizeDM; ++c) {
          std::cout << "    Checking " << faceConeDM[c] << " against " << faceVerticesDM[c] << std::endl;
        } // for
      } // if
      assert(found);
    } // else
    noReplaceCells.insert(otherCell);
    replaceCells.insert(cell);
    if (faceCone) replaceVertices.insert(faceCone, &faceCone[coneSize]);
    replaceVerticesDM.insert(faceConeDM, &faceConeDM[coneSizeDM]);
    cellRenumber[cell]   = firstFaultCell;
    cellRenumberDM[cell] = firstFaultCellDM;
    // Adding cohesive cell (not interpolated)
    PetscInt newv = 0;
    if (debug) std::cout << "  Creating cohesive cell " << firstFaultCell << std::endl;
    for (int c = 0; c < coneSizeDM; ++c) {
      if (debug) std::cout << "    vertex " << faceConeDM[c] << std::endl;
      if (faceCone) sieve->addArrow(faceCone[c], firstFaultCell);
      cohesiveCone[newv++] = faceConeDM[c] + extraCells;
    } // for
    for (int c = 0; c < coneSizeDM; ++c) {
      if (debug) std::cout << "    shadow vertex " << vertexRenumberDM[faceConeDM[c] + extraCells] << std::endl;
      if (faceCone) sieve->addArrow(vertexRenumber[faceCone[c]], firstFaultCell, true);
      cohesiveCone[newv++] = vertexRenumberDM[faceConeDM[c] + extraCells];
    } // for
    if (constraintCell) {
      for (int c = 0; c < coneSizeDM; ++c) {
        if (debug) std::cout << "    Lagrange vertex " << vertexLagrangeRenumberDM[faceConeDM[c] + extraCells] << std::endl;
        if (faceCone) sieve->addArrow(vertexLagrangeRenumber[faceCone[c]], firstFaultCell, true);
        cohesiveCone[newv++] = vertexLagrangeRenumberDM[faceConeDM[c] + extraCells];
      } // for
    } // if
#ifdef USE_DMCOMPLEX_ON
    err = DMPlexSetCone(newMesh, firstFaultCellDM, cohesiveCone);PYLITH_CHECK_ERROR(err);
#endif
    // TODO: Need to reform the material label when sieve is reallocated
    if (faceCone) sieveMesh->setValue(material, firstFaultCell, materialId);
    err = DMPlexSetLabelValue(newMesh, "material-id", firstFaultCellDM, materialId);PYLITH_CHECK_ERROR(err);
    logger.stagePop();
    logger.stagePush("SerialFaultStratification");
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultCell, 0);
    sieveMesh->setDepth(firstFaultCell, 1);
#endif
    logger.stagePop();
    logger.stagePush("SerialFaultCreation");
    if (sV2) {
      sV2->clear();
      cV2->clear();
    }
    err = DMPlexRestoreTransitiveClosure(faultDMMesh, faceDM, PETSC_TRUE, &closureSize, &faceConeDM);PYLITH_CHECK_ERROR(err);
  } // for over fault faces
  delete f_iter; f_iter = NULL;
  if (sV2) {
    delete sV2; sV2 = NULL;
    delete cV2; cV2 = NULL;
  }
  // This completes the set of cells scheduled to be replaced
  // TODO: Convert to DMPlex
  TopologyOps::PointSet replaceCellsBase(replaceCells);

  TopologyOps::PointSet faultBdVertices;
#if 0
  if (!faultBoundary.isNull()) {
    const ALE::Obj<SieveFlexMesh::label_sequence>& faultBdVerts = faultBoundary->depthStratum(0);
    assert(!faultBdVerts.isNull());
    faultBdVertices.insert(faultBdVerts->begin(), faultBdVerts->end());
  } else {
#else
    IS              bdSubpointIS;
    const PetscInt *points;
    PetscInt        bfvStart, bfvEnd;

    assert(faultBoundaryDM);
    err = DMPlexGetDepthStratum(faultBoundaryDM, 0, &bfvStart, &bfvEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexCreateSubpointIS(faultBoundaryDM, &bdSubpointIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(bdSubpointIS, &points);PYLITH_CHECK_ERROR(err);
    for (PetscInt v = bfvStart; v < bfvEnd; ++v) {
      faultBdVertices.insert(subpointMap[points[v]]);
    }
    err = ISRestoreIndices(bdSubpointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&bdSubpointIS);PYLITH_CHECK_ERROR(err);
#endif
  err = ISRestoreIndices(subpointIS, &subpointMap);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);

  TopologyOps::PointSet::const_iterator rVerticesEnd = replaceVerticesDM.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVerticesDM.begin(); v_iter != rVerticesEnd; ++v_iter) {
    if (faultBdVertices.find(*v_iter) != faultBdVertices.end())
      continue;
    //TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, replaceCells, noReplaceCells, debug);
    TopologyOps::classifyCellsDM(complexMesh, *v_iter, depth, faceSizeDM, firstCohesiveCellDM, replaceCells, noReplaceCells, debug);
  } // for
  const TopologyOps::PointSet::const_iterator fbdVerticesEnd = faultBdVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = faultBdVertices.begin(); v_iter != fbdVerticesEnd; ++v_iter) {
    //TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, replaceCells, noReplaceCells, debug);
    TopologyOps::classifyCellsDM(complexMesh, *v_iter, depth, faceSizeDM, firstCohesiveCellDM, replaceCells, noReplaceCells, debug);
  } // for
  
  // Add new arrows for support of replaced vertices
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV(std::max(1, sieve->getMaxSupportSize()));

  rVerticesEnd = replaceVertices.end();
  for(TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != rVerticesEnd; ++v_iter) {
    sieve->support(*v_iter, sV);
    const point_type *support = sV.getPoints();

    if (debug) std::cout << "  Checking support of " << *v_iter << std::endl;
    const int sVSize = sV.getSize();
    for (int s = 0; s < sVSize; ++s) {
      if (replaceCells.find(support[s]) != replaceCells.end()) {
        if (debug) std::cout << "    Adding new support " << vertexRenumber[*v_iter] << " --> " << support[s] << std::endl;
        sieve->addArrow(vertexRenumber[*v_iter], support[s], true);
      } else {
        if (debug) std::cout << "    Keeping same support " << *v_iter<<"," << vertexRenumber[*v_iter] << " --> " << support[s] << std::endl;
      } // if/else
    } // for
    sV.clear();
  }
#ifdef USE_DMCOMPLEX_ON
  TopologyOps::PointSet::const_iterator rVerticesDMEnd = replaceVerticesDM.end();
  for(PetscInt cell = cStart; cell < cEnd; ++cell) {
    const PetscInt *cone;
    PetscInt        coneSize;

    err = DMPlexGetCone(complexMesh, cell, &cone);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(complexMesh, cell, &coneSize);PYLITH_CHECK_ERROR(err);
    if (replaceCells.find(cell) != replaceCells.end()) {
      for(PetscInt c = 0; c < coneSize; ++c) {
        PetscBool replaced = PETSC_FALSE;

        for(TopologyOps::PointSet::const_iterator v_iter = replaceVerticesDM.begin(); v_iter != rVerticesDMEnd; ++v_iter) {
          if (cone[c] == *v_iter) {
            cohesiveCone[c] = vertexRenumberDM[cone[c] + extraCells];
            replaced        = PETSC_TRUE;
            break;
          }
        }
        if (!replaced) {
          cohesiveCone[c] = cone[c] + extraCells;
        }
      }
    } else {
      for(PetscInt c = 0; c < coneSize; ++c) {
        cohesiveCone[c] = cone[c] + extraCells;
      }
    }
    err = DMPlexSetCone(newMesh, cell, cohesiveCone);PYLITH_CHECK_ERROR(err);
  }
#endif
  err = PetscFree3(origVerticesDM, faceVerticesDM, cohesiveCone);PYLITH_CHECK_ERROR(err);
  sieve->reallocate();
#ifdef USE_DMCOMPLEX_ON
  /* DMPlex */
  err = DMPlexSymmetrize(newMesh);PYLITH_CHECK_ERROR(err);
  err = DMPlexStratify(newMesh);PYLITH_CHECK_ERROR(err);
#endif

  // More checking
  const bool firstFault = !sieveMesh->hasRealSection("replaced_cells");
  const ALE::Obj<topology::Mesh::RealSection>& replacedCells = sieveMesh->getRealSection("replaced_cells");
  assert(!replacedCells.isNull());
  TopologyOps::PointSet cellNeighbors;
	 
  if (firstFault) {
    const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->heightStratum(0);
    assert(!cells.isNull());

    replacedCells->setChart(topology::Mesh::RealSection::chart_type(*std::min_element(cells->begin(), cells->end()), *std::max_element(cells->begin(), cells->end())+1));
    replacedCells->setFiberDimension(cells, 1);
    replacedCells->allocatePoint();
  } // if
	 
  const TopologyOps::PointSet::const_iterator noRCellsEnd = noReplaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = noReplaceCells.begin(); c_iter != noRCellsEnd; ++c_iter) {
    const PylithScalar minusOne = -1.0;
    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &minusOne);
    } else {
      const PylithScalar minusTwo = -2.0;
      replacedCells->updatePoint(*c_iter, &minusTwo);
    } // if/else
  } // for

  TopologyOps::PointSet::const_iterator rCellsEnd = replaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = replaceCells.begin(); c_iter != rCellsEnd; ++c_iter) {
    if (replaceCellsBase.find(*c_iter) != replaceCellsBase.end()) {
      const PylithScalar one = 1.0;
      if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
        replacedCells->updatePoint(*c_iter, &one);
      } else {
        const PylithScalar two = 2.0;
        replacedCells->updatePoint(*c_iter, &two);
      } // if/else
      continue;
    } // if
    const PylithScalar ten = 10.0;
    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &ten);
    } else {
      const PylithScalar twenty = 20.0;
      replacedCells->updatePoint(*c_iter, &twenty);
    } // if/else
    // There should be a way to check for boundary elements
    if (mesh->dimension() == 1) {
      if (cellNeighbors.size() > 2) {
        std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
        throw ALE::Exception("Invalid number of neighbors");
      } // if
    } else if (mesh->dimension() == 2) {
      if (numCorners == 3) {
        if (cellNeighbors.size() > 3) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
	} // if
      } else if (numCorners == 4 || numCorners == 9) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } // if/else
    } else if (mesh->dimension() == 3) {
      if (numCorners == 4) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } else if (numCorners == 8) {
        if (cellNeighbors.size() > 6) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } // if/else
    } // if/else
  } // for
  ReplaceVisitor<SieveMesh::sieve_type,std::map<SieveMesh::point_type,SieveMesh::point_type> > rVc(vertexRenumber, std::max(1, sieve->getMaxConeSize()), debug);
  
  rCellsEnd = replaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = replaceCells.begin(); c_iter != rCellsEnd; ++c_iter) {
    sieve->cone(*c_iter, rVc);
    if (rVc.mappedPoint()) {
      if (debug) std::cout << "  Replacing cell " << *c_iter << std::endl;
      sieve->setCone(rVc.getPoints(), *c_iter);
    } // if
    rVc.clear();
  } // for
  ReplaceVisitor<SieveMesh::sieve_type,std::map<SieveMesh::point_type,SieveMesh::point_type> > rVs(cellRenumber, std::max(1, sieve->getMaxSupportSize()), debug);

  rVerticesEnd = replaceVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != rVerticesEnd; ++v_iter) {
    sieve->support(*v_iter, rVs);
    if (rVs.mappedPoint()) {
      if (debug) std::cout << "  Replacing support for " << *v_iter << std::endl;
      sieve->setSupport(*v_iter, rVs.getPoints());
    } else {
      if (debug) std::cout << "  Not replacing support for " << *v_iter << std::endl;
    } // if/else
    rVs.clear();
  } // for
  if (!rank) {
    delete [] indices;
    err = PetscFree(indicesDM);PYLITH_CHECK_ERROR(err);
  }
#if !defined(FAST_STRATIFY)
  logger.stagePop();
  logger.stagePush("SerialFaultStratification");
  sieveMesh->stratify();
  logger.stagePop();
  logger.stagePush("SerialFaultCreation");
#endif
  const std::string labelName("censored depth");

  if (!sieveMesh->hasLabel(labelName)) {
    const ALE::Obj<SieveMesh::label_type>& label = sieveMesh->createLabel(labelName);
    assert(!label.isNull());

    TopologyOps::computeCensoredDepth(label, sieveMesh->getSieve(), firstFaultVertex);
  } else {
    // Insert new shadow vertices into existing label
    const ALE::Obj<SieveMesh::label_type>& label = sieveMesh->getLabel(labelName);
    assert(!label.isNull());

    const std::map<int,int>::const_iterator vRenumberEnd = vertexRenumber.end();
    for (std::map<int,int>::const_iterator v_iter = vertexRenumber.begin(); v_iter != vRenumberEnd; ++v_iter)
      sieveMesh->setValue(label, v_iter->second, 0);
  } // if/else
  if (debug) mesh->view("Mesh with Cohesive Elements");

  // Fix coordinates
  const ALE::Obj<topology::Mesh::RealSection>& coordinates = sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  PetscSection coordSection, newCoordSection;
  Vec          coordinatesVec, newCoordinatesVec;
  PetscScalar *coords, *newCoords;
  PetscInt     numComp, coordSize;
 
  err = DMPlexGetCoordinateSection(complexMesh, &coordSection);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetCoordinateSection(newMesh,     &newCoordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetNumFields(newCoordSection, 1);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetDof(coordSection, vStart, &numComp);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetFieldComponents(newCoordSection, 0, numComp);PYLITH_CHECK_ERROR(err);
  err = DMGetCoordinatesLocal(complexMesh, &coordinatesVec);PYLITH_CHECK_ERROR(err);
  err = PetscSectionSetChart(newCoordSection, vStart+extraCells, vEnd+extraCells+extraVertices);PYLITH_CHECK_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof;
    err = PetscSectionGetDof(coordSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetDof(newCoordSection, v+extraCells, dof);PYLITH_CHECK_ERROR(err);
  }

  if (debug) coordinates->view("Coordinates without shadow vertices");
  for (PetscInt fv = 0; fv < numFaultVerticesDM; ++fv) {
    PetscInt v    = fVerticesDM[fv];
    PetscInt vnew = v+extraCells;
    PetscInt dof;

    err = PetscSectionGetDof(coordSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetDof(newCoordSection, vertexRenumberDM[vnew], dof);PYLITH_CHECK_ERROR(err);
    if (constraintCell) {err = PetscSectionSetDof(newCoordSection, vertexLagrangeRenumberDM[vnew], dof);PYLITH_CHECK_ERROR(err);}

    coordinates->addPoint(vertexRenumber[v], coordinates->getFiberDimension(v));
    if (constraintCell) coordinates->addPoint(vertexLagrangeRenumber[v], coordinates->getFiberDimension(v));
  } // for
  sieveMesh->reallocate(coordinates);

  err = PetscSectionSetUp(newCoordSection);PYLITH_CHECK_ERROR(err);
  err = PetscSectionGetStorageSize(newCoordSection, &coordSize);PYLITH_CHECK_ERROR(err);
  err = VecCreate(mesh->comm(), &newCoordinatesVec);PYLITH_CHECK_ERROR(err);
  err = VecSetSizes(newCoordinatesVec, coordSize, PETSC_DETERMINE);PYLITH_CHECK_ERROR(err);
  err = VecSetFromOptions(newCoordinatesVec);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(coordinatesVec, &coords);PYLITH_CHECK_ERROR(err);
  err = VecGetArray(newCoordinatesVec, &newCoords);PYLITH_CHECK_ERROR(err);

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt vnew = v+extraCells;
    PetscInt dof, off, newoff, d;

    err = PetscSectionGetDof(coordSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(newCoordSection, vnew, &newoff);PYLITH_CHECK_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      newCoords[newoff+d] = coords[off+d];
    }
  }
  for (PetscInt fv = 0; fv < numFaultVerticesDM; ++fv) {
    PetscInt v    = fVerticesDM[fv];

    assert(coordinates->getFiberDimension(v) == coordinates->getFiberDimension(vertexRenumber[v]));
    coordinates->updatePoint(vertexRenumber[v], coordinates->restrictPoint(v));
    if (constraintCell) {
      assert(coordinates->getFiberDimension(v) == coordinates->getFiberDimension(vertexLagrangeRenumber[v]));
      coordinates->updatePoint(vertexLagrangeRenumber[v], coordinates->restrictPoint(v));
    } // if
  } // for
  for (PetscInt fv = 0; fv < numFaultVerticesDM; ++fv) {
    PetscInt v    = fVerticesDM[fv];
    PetscInt vnew = v+extraCells;
    PetscInt dof, off, newoff, d;

    err = PetscSectionGetDof(coordSection, v, &dof);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetOffset(newCoordSection, vertexRenumberDM[vnew], &newoff);PYLITH_CHECK_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      newCoords[newoff+d] = coords[off+d];
    }
    if (constraintCell) {
      err = PetscSectionGetOffset(newCoordSection, vertexLagrangeRenumberDM[vnew], &newoff);PYLITH_CHECK_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        newCoords[newoff+d] = coords[off+d];
      }
    } // if
  } // for
  err = VecRestoreArray(coordinatesVec, &coords);PYLITH_CHECK_ERROR(err);
  err = VecRestoreArray(newCoordinatesVec, &newCoords);PYLITH_CHECK_ERROR(err);
  err = DMSetCoordinatesLocal(newMesh, newCoordinatesVec);PYLITH_CHECK_ERROR(err);
  err = VecDestroy(&newCoordinatesVec);PYLITH_CHECK_ERROR(err);
  if (debug) coordinates->view("Coordinates with shadow vertices");

  err = ISRestoreIndices(fVertexIS, &fVerticesDM);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&fVertexIS);PYLITH_CHECK_ERROR(err);
  logger.stagePop();

  mesh->setDMMesh(newMesh);

  PYLITH_METHOD_END;
} // create

void
pylith::faults::CohesiveTopology::createInterpolated(topology::Mesh* mesh,
                                         const topology::SubMesh& faultMesh,
                                         const ALE::Obj<SieveFlexMesh>& faultBoundary,
                                         const ALE::Obj<topology::Mesh::IntSection>& groupField,
                                         const int materialId,
                                         int& firstFaultVertex,
                                         int& firstLagrangeVertex,
                                         int& firstFaultCell,
                                         const bool constraintCell)
{ // createInterpolated
  assert(0 != mesh);
  assert(!faultBoundary.isNull());
  assert(!groupField.isNull());
  DM             dm = mesh->dmMesh(), sdm;
  DMLabel        label;
  const char    *labelName = "faultSurface";
  PetscErrorCode err;

  err = DMPlexGetLabel(dm, labelName, &label);PYLITH_CHECK_ERROR(err);
  // Completes the set of cells scheduled to be replaced
  //   Have to do internal fault vertices before fault boundary vertices, and this is the only thing I use faultBoundary for
  err = DMPlexLabelCohesiveComplete(dm, label);PYLITH_CHECK_ERROR(err);
  err = DMPlexConstructCohesiveCells(dm, label, &sdm);PYLITH_CHECK_ERROR(err);
  mesh->setDMMesh(sdm);
} // createInterpolated

PetscInt convertSieveToDMPointNumbering(PetscInt sievePoint, PetscInt numNormalCells, PetscInt numCohesiveCells, PetscInt numNormalVertices, PetscInt numShadowVertices, PetscInt numLagrangeVertices)
{
  PetscInt dmPoint = -1;

  if ((sievePoint >= 0) && (sievePoint < numNormalCells)) {
    dmPoint = sievePoint;
    //std::cout << "normal cell sieve point "<<sievePoint<<" --> "<<" dm point"<<dmPoint<<std::endl;
  } else if ((sievePoint >= numNormalCells) && (sievePoint < numNormalCells+numNormalVertices)) {
    dmPoint = sievePoint+numCohesiveCells;
    //std::cout << "normal vertex sieve point "<<sievePoint<<" --> "<<" dm point"<<dmPoint<<std::endl;
  } else if ((sievePoint >= numNormalCells+numNormalVertices) && (sievePoint < numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices)) {
    dmPoint = sievePoint+numCohesiveCells;
    //std::cout << "extra vertex sieve point "<<sievePoint<<" --> "<<" dm point"<<dmPoint<<std::endl;
  } else if ((sievePoint >= numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices) && (sievePoint < numNormalCells+numNormalVertices+numShadowVertices+numLagrangeVertices+numCohesiveCells)) {
    dmPoint = sievePoint-(numNormalVertices+numShadowVertices+numLagrangeVertices);
    //std::cout << "extra cell sieve point "<<sievePoint<<" --> "<<" dm point"<<dmPoint<<std::endl;
  } else {
    //std::cout << "face sieve point "<<sievePoint<<" --> "<<" dm point"<<dmPoint<<std::endl;
  }
  return dmPoint;
}

PetscInt convertDMToSievePointNumbering(PetscInt dmPoint, PetscInt numNormalCells, PetscInt numCohesiveCells, PetscInt numNormalVertices, PetscInt numShadowVertices, PetscInt numLagrangeVertices)
{
  PetscInt sievePoint = -1;

  if ((dmPoint >= 0) && (dmPoint < numNormalCells)) {
    sievePoint = dmPoint;
    //std::cout << "normal cell sieve point "<<sievePoint<<" <-- "<<" dm point"<<dmPoint<<std::endl;
  } else if ((dmPoint >= numNormalCells) && (dmPoint < numNormalCells+numCohesiveCells)) {
    sievePoint = dmPoint+numNormalVertices+numShadowVertices+numLagrangeVertices;
    //std::cout << "extra cell sieve point "<<sievePoint<<" <-- "<<" dm point"<<dmPoint<<std::endl;
  } else if ((dmPoint >= numNormalCells+numCohesiveCells) && (dmPoint < numNormalCells+numCohesiveCells+numNormalVertices)) {
    sievePoint = dmPoint-numCohesiveCells;
    //std::cout << "normal vertex sieve point "<<sievePoint<<" <-- "<<" dm point"<<dmPoint<<std::endl;
  } else if ((dmPoint >= numNormalCells+numCohesiveCells+numNormalVertices) && (dmPoint < numNormalCells+numCohesiveCells+numNormalVertices+numShadowVertices+numLagrangeVertices)) {
    sievePoint = dmPoint-numCohesiveCells;
    //std::cout << "extra vertex sieve point "<<sievePoint<<" <-- "<<" dm point"<<dmPoint<<std::endl;
  } else {
    //std::cout << "face sieve point "<<sievePoint<<" <-- "<<" dm point"<<dmPoint<<std::endl;
  }
  return sievePoint;
}

// ----------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::CohesiveTopology::createFaultParallel(
			    topology::SubMesh* faultMesh,
			    const topology::Mesh& mesh,
			    const int materialId,
			    const bool constraintCell)
{ // createFaultParallel
  PYLITH_METHOD_BEGIN;

  assert(faultMesh);

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultCreation");
  PetscErrorCode err;

  faultMesh->coordsys(mesh);

  DM dmMesh = mesh.dmMesh();
  assert(dmMesh);
  DM dmFaultMesh;

  err = DMPlexCreateCohesiveSubmesh(dmMesh, constraintCell ? PETSC_TRUE : PETSC_FALSE, &dmFaultMesh);PYLITH_CHECK_ERROR(err);
  faultMesh->setDMMesh(dmFaultMesh);

  logger.stagePop();

  PYLITH_METHOD_END;
} // createFaultParallel


// End of file
