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

#include "CohesiveTopology.hh" // implementation of object methods

#include "TopologyOps.hh" // USES TopologyOps
#include "TopologyVisitors.hh" // USES TopologyVisitors
#include "pylith/topology/SubMesh.hh" // USES SubMesh

#include <Selection.hh> // Algorithms for submeshes

#include <cassert> // USES assert()

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

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
                                              const topology::Mesh& mesh,
                                              const ALE::Obj<topology::Mesh::IntSection>& groupField,
                                              const bool flipFault)
{ // createFault
  assert(0 != faultMesh);
  assert(!groupField.isNull());

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("SerialFaultCreation");
  PetscErrorCode err;

  faultMesh->coordsys(mesh);
  DM       dmMesh = mesh.dmMesh();
  PetscInt dim, depth;

  assert(dmMesh);
  err = DMPlexGetDimension(dmMesh, &dim);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepth(dmMesh, &depth);CHECK_PETSC_ERROR(err);

  const ALE::Obj<SieveMesh>&             sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::sieve_type>& sieve     = sieveMesh->getSieve();
  assert(!sieve.isNull());
  ALE::Obj<SieveSubMesh>&                faultSieveMesh = faultMesh->sieveMesh();
  faultSieveMesh = new SieveSubMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());
  const ALE::Obj<SieveSubMesh::sieve_type> ifaultSieve = new SieveMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!ifaultSieve.isNull());

  ALE::Obj<SieveFlexMesh>             fault      = new SieveFlexMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());
  ALE::Obj<SieveFlexMesh::sieve_type> faultSieve = new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!fault.isNull());assert(!faultSieve.isNull());
  const int debug = mesh.debug();

  // Create set with vertices on fault
  const IntSection::chart_type& chart = groupField->getChart();
  TopologyOps::PointSet faultVertices; // Vertices on fault

  const IntSection::chart_type::const_iterator chartEnd = chart.end();
  for(IntSection::chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    assert(!sieveMesh->depth(*c_iter));
    if (groupField->getFiberDimension(*c_iter))
      faultVertices.insert(*c_iter);
  } // for

  // Create a sieve which captures the fault
  const bool vertexFault = true;
  const int firstFaultCell = sieve->getBaseSize() + sieve->getCapSize();

  TopologyOps::createFaultSieveFromVertices(fault->getDimension(), firstFaultCell, 
					    faultVertices, sieveMesh, 
					    fault->getArrowSection("orientation"), 
					    faultSieve, flipFault);
  fault->setSieve(faultSieve);

  logger.stagePop();
  logger.stagePush("SerialFaultStratification");
  fault->stratify();
  logger.stagePop();
  logger.stagePush("SerialFaultCreation");
  if (debug)
    fault->view("Fault mesh");

  faultBoundary = ALE::Selection<SieveFlexMesh>::boundary(fault);
  if (debug)
    faultBoundary->view("Fault boundary mesh");

  // Orient the fault sieve
  TopologyOps::orientFaultSieve(fault->getDimension(), sieveMesh,
				fault->getArrowSection("orientation"), fault);

  // Convert fault to an IMesh
  SieveSubMesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
  faultSieveMesh->setSieve(ifaultSieve);
  ALE::ISieveConverter::convertMesh(*fault, *faultSieveMesh, renumbering, false);
  renumbering.clear();

  // Convert fault to a DM
  if (depth == dim) {
    DM             subdm;
    DMLabel        label;
    const char    *labelName = groupField->getName().c_str();

    // Put fault vertices in a label
    err = DMPlexCreateLabel(dmMesh, labelName);CHECK_PETSC_ERROR(err);
    err = DMPlexGetLabel(dmMesh, labelName, &label);CHECK_PETSC_ERROR(err);
    const IntSection::chart_type& chart = groupField->getChart();
    const IntSection::chart_type::const_iterator chartEnd = chart.end();
    for(IntSection::chart_type::const_iterator c_iter = chart.begin(); c_iter != chartEnd; ++c_iter) {
      if (groupField->getFiberDimension(*c_iter)) {
        err = DMLabelSetValue(label, *c_iter, 0);CHECK_PETSC_ERROR(err);
      }
    }

    err = DMPlexCreateSubmesh(dmMesh, labelName, &subdm);CHECK_PETSC_ERROR(err);
    faultMesh->setDMMesh(subdm);
  } else {
    // TODO: This leg will be unnecessary
    DM             dm;
    DMLabel        subpointMap;
    PetscInt      *renum;
    PetscInt       pStart, pEnd, cStart, cEnd, vStart, vEnd;
    PetscErrorCode err;
    SieveSubMesh::renumbering_type renumbering;

    ALE::ISieveConverter::convertMesh(*fault, &dm, renumbering, true);
    // Have to make subpointMap here: renumbering[original] = fault
    err = DMLabelCreate("subpoint_map", &subpointMap);CHECK_PETSC_ERROR(err);
    err = DMPlexGetChart(dm, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    assert(renumbering.size() == pEnd-pStart);
    err = PetscMalloc((pEnd-pStart) * sizeof(PetscInt), &renum);CHECK_PETSC_ERROR(err);
    for(SieveSubMesh::renumbering_type::const_iterator p_iter = renumbering.begin(); p_iter != renumbering.end(); ++p_iter) {
      renum[p_iter->second] = p_iter->first;
#if 0
      std::cout << "renum["<<p_iter->second<<"]: "<<p_iter->first<<std::endl;
#endif
    }
    for(PetscInt p = 1; p < pEnd-pStart; ++p) {
      assert(renum[p] > renum[p-1]);
    }
    err = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    for(PetscInt p = cStart; p < cEnd; ++p) {
      err = DMLabelSetValue(subpointMap, renum[p], mesh.dimension());CHECK_PETSC_ERROR(err);
    }
    err = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
    for(PetscInt p = vStart; p < vEnd; ++p) {
      err = DMLabelSetValue(subpointMap, renum[p], 0);CHECK_PETSC_ERROR(err);
    }
    err = PetscFree(renum);CHECK_PETSC_ERROR(err);
    err = DMPlexSetSubpointMap(dm, subpointMap);CHECK_PETSC_ERROR(err);
    renumbering.clear();
    faultMesh->setDMMesh(dm);
  }

  logger.stagePop();
} // createFault

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::create(topology::Mesh* mesh,
                                         const topology::SubMesh& faultMesh,
                                         const ALE::Obj<SieveFlexMesh>& faultBoundary,
                                         const ALE::Obj<topology::Mesh::IntSection>& groupField,
                                         const int materialId,
                                         int& firstFaultVertex,
                                         int& firstLagrangeVertex,
                                         int& firstFaultCell,
                                         const bool constraintCell)
{ // create
  assert(0 != mesh);
  assert(!faultBoundary.isNull());
  assert(!groupField.isNull());

  typedef ALE::SieveAlg<SieveFlexMesh> sieveAlg;
  typedef ALE::Selection<SieveFlexMesh> selection;
  PetscErrorCode err;

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("SerialFaultCreation");

  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  /* DMPlex */
  DM complexMesh = mesh->dmMesh();
  assert(complexMesh);
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  assert(!faultSieveMesh.isNull());  

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<SieveSubMesh::sieve_type> ifaultSieve = faultSieveMesh->getSieve();
  assert(!ifaultSieve.isNull());

  const int  depth = sieveMesh->depth();
  assert(!sieveMesh->heightStratum(0).isNull());
  const int numCells = sieveMesh->heightStratum(0)->size();
  PetscInt cStart, cEnd;
#define USE_DMCOMPLEX_ON
#ifdef USE_DMCOMPLEX_ON
  err = DMPlexGetHeightStratum(complexMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  assert(numCells == cEnd - cStart);
#endif
  int numCorners = 0; // The number of vertices in a mesh cell
  int faceSize = 0; // The number of vertices in a mesh face
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

  if (!faultSieveMesh->commRank()) {
    assert(!faultSieveMesh->heightStratum(1).isNull());
    const SieveSubMesh::point_type p = *faultSieveMesh->heightStratum(1)->begin();

    numCorners = sieveMesh->getNumCellCorners();
    faceSize = selection::numFaceVertices(sieveMesh);
    indices = new int[faceSize];
    numFaultCorners = faultSieveMesh->getNumCellCorners(p, faultSieveMesh->depth(p));

#ifdef USE_DMCOMPLEX_ON
    PetscInt numCornersDM;
    err = DMPlexGetConeSize(complexMesh, cStart, &numCornersDM);CHECK_PETSC_ERROR(err);
    assert(numCorners == numCornersDM);
    err = PetscMalloc(faceSize * sizeof(PetscInt), &indicesDM);CHECK_PETSC_ERROR(err);
#endif
    /* TODO: Do we need faceSize at all? Blaise was using a nice criterion */
  }
  //faultSieveMesh->view("Serial fault mesh");

  // Add new shadow vertices and possibly Lagrange multipler vertices
  const ALE::Obj<SieveSubMesh::label_sequence>& fVertices           = faultSieveMesh->depthStratum(0);
  assert(!fVertices.isNull());
  const SieveSubMesh::label_sequence::const_iterator fVerticesBegin = fVertices->begin();
  const SieveSubMesh::label_sequence::const_iterator fVerticesEnd   = fVertices->end();
  const ALE::Obj<SieveMesh::label_sequence>&         vertices       = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
#ifdef USE_DMCOMPLEX_ON
  PetscInt vStart, vEnd;
  err = DMPlexGetDepthStratum(complexMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  assert(vertices->size() == vEnd - vStart);
#endif
  const ALE::Obj<std::set<std::string> >& groupNames = sieveMesh->getIntSections();
  assert(!groupNames.isNull());
  const int numFaultVertices = fVertices->size();
  std::map<point_type,point_type> vertexRenumber;
  std::map<point_type,point_type> vertexLagrangeRenumber;
  std::map<point_type,point_type> cellRenumber;
  std::map<point_type,point_type> vertexRenumberDM;
  std::map<point_type,point_type> vertexLagrangeRenumberDM;
  std::map<point_type,point_type> cellRenumberDM;
  PetscInt numFaultFaces = faultSieveMesh->heightStratum(1)->size();
  PetscInt faultVertexOffsetDM, firstFaultVertexDM, firstLagrangeVertexDM, firstFaultCellDM, extraVertices, extraCells;
#ifdef USE_DMCOMPLEX_ON
  /* TODO This will have to change for multiple faults */
  PetscInt numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices;

  extraVertices         = numFaultVertices * (constraintCell ? 2 : 1); /* Total number of fault vertices on this fault (shadow + Lagrange) */
  extraCells            = numFaultFaces;                               /* Total number of fault cells */
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
    err = DMPlexGetChart(complexMesh, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    assert(firstFaultVertex == pEnd - pStart);
#endif
  }

#ifdef USE_DMCOMPLEX_ON
  /* DMPlex */
  DM        newMesh;
  PetscInt *newCone;
  PetscInt  maxConeSize = 0;

  err = DMCreate(sieveMesh->comm(), &newMesh);CHECK_PETSC_ERROR(err);
  err = DMSetType(newMesh, DMPLEX);CHECK_PETSC_ERROR(err);
  err = DMPlexSetDimension(newMesh, sieveMesh->getDimension());CHECK_PETSC_ERROR(err);
  err = DMPlexSetChart(newMesh, 0, firstFaultVertexDM + extraVertices);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt coneSize;
    err = DMPlexGetConeSize(complexMesh, c, &coneSize);CHECK_PETSC_ERROR(err);
    err = DMPlexSetConeSize(newMesh, c, coneSize);CHECK_PETSC_ERROR(err);
    maxConeSize = PetscMax(maxConeSize, coneSize);
  }
  for(PetscInt c = cEnd; c < cEnd+numFaultFaces; ++c) {
    err = DMPlexSetConeSize(newMesh, c, constraintCell ? faceSize*3 : faceSize*2);CHECK_PETSC_ERROR(err);
  }
  err = DMSetUp(newMesh);CHECK_PETSC_ERROR(err);
  err = PetscMalloc(maxConeSize * sizeof(PetscInt), &newCone);CHECK_PETSC_ERROR(err);
  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt *cone;
    PetscInt        coneSize, cp;

    err = DMPlexGetCone(complexMesh, c, &cone);CHECK_PETSC_ERROR(err);
    err = DMPlexGetConeSize(complexMesh, c, &coneSize);CHECK_PETSC_ERROR(err);
    for(cp = 0; cp < coneSize; ++cp) {
      newCone[cp] = cone[cp] + extraCells;
    }
    err = DMPlexSetCone(newMesh, c, newCone);CHECK_PETSC_ERROR(err);
  }
  err = PetscFree(newCone);CHECK_PETSC_ERROR(err);
  PetscInt cMax, vMax;

  err = DMPlexGetHybridBounds(newMesh, &cMax, PETSC_NULL, PETSC_NULL, &vMax);CHECK_PETSC_ERROR(err);
  if (cMax < 0) {
    err = DMPlexSetHybridBounds(newMesh, firstFaultCellDM, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  }
  err = DMPlexSetHybridBounds(newMesh, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DETERMINE, firstLagrangeVertexDM);CHECK_PETSC_ERROR(err);

  // Renumber labels
  std::set<std::string> names(groupNames->begin(), groupNames->end());
  names.insert(names.begin(), "material-id");
  for(std::set<std::string>::const_iterator name = names.begin(); name != names.end(); ++name) {
    const char     *lname = (*name).c_str();
    IS              idIS;
    PetscInt        n;
    PetscBool       hasLabel;
    const PetscInt *ids;

    err = DMPlexHasLabel(complexMesh, lname, &hasLabel);CHECK_PETSC_ERROR(err);
    if (!hasLabel) continue;
    err = DMPlexGetLabelSize(complexMesh, lname, &n);CHECK_PETSC_ERROR(err);
    err = DMPlexGetLabelIdIS(complexMesh, lname, &idIS);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(idIS, &ids);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < n; ++i) {
      const PetscInt  id = ids[i];
      const PetscInt *points;
      IS              sIS;
      PetscInt        size;

      err = DMPlexGetStratumSize(complexMesh, lname, id, &size);CHECK_PETSC_ERROR(err);
      err = DMPlexGetStratumIS(complexMesh, lname, id, &sIS);CHECK_PETSC_ERROR(err);
      err = ISGetIndices(sIS, &points);CHECK_PETSC_ERROR(err);
      for(PetscInt s = 0; s < size; ++s) {
        if ((points[s] >= vStart) && (points[s] < vEnd)) {
          err = DMPlexSetLabelValue(newMesh, lname, points[s]+extraCells, id);CHECK_PETSC_ERROR(err);
        } else {
          err = DMPlexSetLabelValue(newMesh, lname, points[s], id);CHECK_PETSC_ERROR(err);
        }
      }
      err = ISRestoreIndices(sIS, &points);CHECK_PETSC_ERROR(err);
      err = ISDestroy(&sIS);CHECK_PETSC_ERROR(err);
    }
    err = ISRestoreIndices(idIS, &ids);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&idIS);CHECK_PETSC_ERROR(err);
  } // for
#endif

  // Add fault vertices to groups and construct renumberings
  for(SieveSubMesh::label_sequence::iterator v_iter = fVerticesBegin; v_iter != fVerticesEnd; ++v_iter, ++firstFaultVertex, ++firstFaultVertexDM) {
    vertexRenumber[*v_iter] = firstFaultVertex;
    vertexRenumberDM[*v_iter+faultVertexOffsetDM] = firstFaultVertexDM;
    if (debug) std::cout << "Duplicating " << *v_iter << " to "	<< vertexRenumber[*v_iter] << std::endl;

    logger.stagePop();
    logger.stagePush("SerialFaultStratification");
    // Add shadow and constraint vertices (if they exist) to group
    // associated with fault
    groupField->addPoint(firstFaultVertex, 1);
    err = DMPlexSetLabelValue(newMesh, groupField->getName().c_str(), firstFaultVertexDM, 1);CHECK_PETSC_ERROR(err);
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultVertex, 1);
    sieveMesh->setDepth(firstFaultVertex, 0);
#endif
    if (constraintCell) {
      vertexLagrangeRenumber[*v_iter] = firstLagrangeVertex;
      vertexLagrangeRenumberDM[*v_iter+faultVertexOffsetDM] = firstLagrangeVertexDM;
      groupField->addPoint(firstLagrangeVertex, 1);
      err = DMPlexSetLabelValue(newMesh, groupField->getName().c_str(), firstLagrangeVertexDM, 1);CHECK_PETSC_ERROR(err);
#if defined(FAST_STRATIFY)
      // OPTIMIZATION
      sieveMesh->setHeight(firstLagrangeVertex, 1);
      sieveMesh->setDepth(firstLagrangeVertex, 0);
#endif
      ++firstLagrangeVertex;
      ++firstLagrangeVertexDM;
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
      if (group->getFiberDimension(*v_iter))
        group->addPoint(firstFaultVertex, 1);

      PetscInt vertexDM = *v_iter+faultVertexOffsetDM;
      PetscInt value;
      //assert(extraCells == faultVertexOffsetDM);
      err = DMPlexGetLabelValue(complexMesh, (*name).c_str(), vertexDM, &value);CHECK_PETSC_ERROR(err);
      if (value) {
        err = DMPlexSetLabelValue(newMesh, (*name).c_str(), vertexRenumberDM[vertexDM], value);CHECK_PETSC_ERROR(err);
      }
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
  const ALE::Obj<SieveSubMesh::label_sequence>& faces = faultSieveMesh->heightStratum(1);
  assert(!faces.isNull());
  const SieveSubMesh::label_sequence::const_iterator facesBegin = faces->begin();
  const SieveSubMesh::label_sequence::const_iterator facesEnd = faces->end();
  const ALE::Obj<SieveFlexMesh::label_type>& material = sieveMesh->getLabel("material-id");
  assert(!material.isNull());
  const int firstCohesiveCell = firstFaultCell;
  TopologyOps::PointSet replaceCells;
  TopologyOps::PointSet noReplaceCells;
  TopologyOps::PointSet replaceVertices;
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV2(std::max(1, ifaultSieve->getMaxSupportSize()));
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> cV2(*ifaultSieve, (size_t) pow(std::max(1, ifaultSieve->getMaxConeSize()), faultSieveMesh->depth()));
  std::set<SieveFlexMesh::point_type> faceSet;
  PetscInt *cohesiveCone;

  err = PetscMalloc3(faceSize,PetscInt,&origVerticesDM,faceSize,PetscInt,&faceVerticesDM,faceSize*3,PetscInt,&cohesiveCone);CHECK_PETSC_ERROR(err);
  for(SieveSubMesh::label_sequence::iterator f_iter = facesBegin; f_iter != facesEnd; ++f_iter, ++firstFaultCell, ++firstFaultCellDM) {
    const point_type face = *f_iter;
    const PetscInt faceConeSizeDM = 10;
    PetscInt faceConeDM[10];
    if (debug) std::cout << "Considering fault face " << face << std::endl;
    ifaultSieve->support(face, sV2);
    const point_type *cells = sV2.getPoints();
    point_type cell = cells[0];
    point_type otherCell = cells[1];

    if (debug) std::cout << "  Checking orientation against cell " << cell << std::endl;
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*ifaultSieve, face, cV2);
    const int coneSize = cV2.getSize();
    const point_type *faceCone = cV2.getPoints();
    //ifaultSieve->cone(face, cV2);
    //const int coneSize = cV2.getSize() ? cV2.getSize() : 1;
    //const point_type *faceCone = cV2.getSize() ? cV2.getPoints() : &face;
    bool found = true;

    for(int i = 0; i < coneSize; ++i) {
      faceSet.insert(faceCone[i]);
      faceConeDM[i] = faceCone[i]+faultVertexOffsetDM;
    }
    selection::getOrientedFace(sieveMesh, cell, &faceSet, numCorners, indices, &origVertices, &faceVertices);
    if (faceVertices.size() != coneSize) {
      std::cout << "Invalid size for faceVertices " << faceVertices.size() << " of face " << face << "should be " << coneSize << std::endl;
      std::cout << "  firstCohesiveCell " << firstCohesiveCell << " firstFaultCell " << firstFaultCell << " numFaces " << faces->size() << std::endl;
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
    assert(faceConeSizeDM >= coneSize);
    assert(faceVertices.size() == coneSize);
    faceSet.clear();
#ifdef USE_DMCOMPLEX_ON
    err = DMPlexGetOrientedFace(complexMesh, cell, coneSize, faceConeDM, numCorners, indicesDM, origVerticesDM, faceVerticesDM, PETSC_NULL);CHECK_PETSC_ERROR(err);
    for(PetscInt c = 0; c < coneSize; ++c) {
      assert(faceVertices[c]+faultVertexOffsetDM == faceVerticesDM[c]);
    }
#endif

    if (numFaultCorners == 0) {
      found = false;
    } else if (numFaultCorners == 2) {
      if (faceVertices[0] != faceCone[0])
        found = false;
    } else {
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVertices[v] != faceCone[0]))
        ++v;
      for(int c = 0; c < coneSize; ++c, ++v) {
        if (debug) std::cout << "    Checking " << faceCone[c] << " against " << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != faceCone[c]) {
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
        while((v < numFaultCorners) && (faceVertices[v] != faceCone[coneSize-1]))
          ++v;
        for(int c = coneSize-1; c >= 0; --c, ++v) {
          if (debug) std::cout << "    Checking " << faceCone[c] << " against " << faceVertices[v%numFaultCorners] << std::endl;
          if (faceVertices[v%numFaultCorners] != faceCone[c]) {
            found = false;
            break;
          } // if
        } // for
      } // if
      if (!found) {
        std::cout << "Considering fault face " << face << std::endl;
        std::cout << "  bordered by cells " << cell << " and " << otherCell << std::endl;
        for(int c = 0; c < coneSize; ++c) {
          std::cout << "    Checking " << faceCone[c] << " against " << faceVertices[c] << std::endl;
        } // for
      } // if
      assert(found);
    } // else
    noReplaceCells.insert(otherCell);
    replaceCells.insert(cell);
    replaceVertices.insert(faceCone, &faceCone[coneSize]);
    cellRenumber[cell]   = firstFaultCell;
    cellRenumberDM[cell] = firstFaultCellDM;
    // Adding cohesive cell (not interpolated)
    PetscInt newv = 0;
    if (debug) std::cout << "  Creating cohesive cell " << firstFaultCell << std::endl;
    for (int c = 0; c < coneSize; ++c) {
      if (debug) std::cout << "    vertex " << faceCone[c] << std::endl;
      sieve->addArrow(faceCone[c], firstFaultCell);
      cohesiveCone[newv++] = faceConeDM[c] + extraCells;
    } // for
    for (int c = 0; c < coneSize; ++c) {
      if (debug) std::cout << "    shadow vertex " << vertexRenumber[faceCone[c]] << std::endl;
      sieve->addArrow(vertexRenumber[faceCone[c]], firstFaultCell, true);
      cohesiveCone[newv++] = vertexRenumberDM[faceConeDM[c]];
    } // for
    if (constraintCell) {
      for (int c = 0; c < coneSize; ++c) {
        if (debug) std::cout << "    Lagrange vertex " << vertexLagrangeRenumber[faceCone[c]] << std::endl;
        sieve->addArrow(vertexLagrangeRenumber[faceCone[c]], firstFaultCell, true);
        cohesiveCone[newv++] = vertexLagrangeRenumberDM[faceConeDM[c]];
      } // for
    } // if
#ifdef USE_DMCOMPLEX_ON
    err = DMPlexSetCone(newMesh, firstFaultCellDM, cohesiveCone);CHECK_PETSC_ERROR(err);
#endif
    // TODO: Need to reform the material label when sieve is reallocated
    sieveMesh->setValue(material, firstFaultCell, materialId);
    err = DMPlexSetLabelValue(newMesh, "material-id", firstFaultCellDM, materialId);CHECK_PETSC_ERROR(err);
    logger.stagePop();
    logger.stagePush("SerialFaultStratification");
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultCell, 0);
    sieveMesh->setDepth(firstFaultCell, 1);
#endif
    logger.stagePop();
    logger.stagePush("SerialFaultCreation");
    sV2.clear();
    cV2.clear();
  } // for over fault faces
  // This completes the set of cells scheduled to be replaced
  // TODO: Convert to DMPlex
  TopologyOps::PointSet replaceCellsBase(replaceCells);

  const ALE::Obj<SieveFlexMesh::label_sequence>& faultBdVerts = faultBoundary->depthStratum(0);
  assert(!faultBdVerts.isNull());
  TopologyOps::PointSet faultBdVertices;

  faultBdVertices.insert(faultBdVerts->begin(), faultBdVerts->end());
  TopologyOps::PointSet::const_iterator rVerticesEnd = replaceVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != rVerticesEnd; ++v_iter) {
    if (faultBdVertices.find(*v_iter) != faultBdVertices.end())
      continue;
    TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, replaceCells, noReplaceCells, debug);
  } // for
  const TopologyOps::PointSet::const_iterator fbdVerticesEnd = faultBdVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = faultBdVertices.begin(); v_iter != fbdVerticesEnd; ++v_iter) {
    TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, replaceCells, noReplaceCells, debug);
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
  for(PetscInt cell = cStart; cell < cEnd; ++cell) {
    const PetscInt *cone;
    PetscInt        coneSize;

    err = DMPlexGetCone(complexMesh, cell, &cone);CHECK_PETSC_ERROR(err);
    err = DMPlexGetConeSize(complexMesh, cell, &coneSize);CHECK_PETSC_ERROR(err);
    if (replaceCells.find(cell) != replaceCells.end()) {
      for(PetscInt c = 0; c < coneSize; ++c) {
        PetscBool replaced = PETSC_FALSE;

        for(TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != rVerticesEnd; ++v_iter) {
          if (cone[c] == *v_iter) {
            cohesiveCone[c] = vertexRenumberDM[cone[c]];
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
    err = DMPlexSetCone(newMesh, cell, cohesiveCone);CHECK_PETSC_ERROR(err);
  }
#endif
  err = PetscFree3(origVerticesDM, faceVerticesDM, cohesiveCone);CHECK_PETSC_ERROR(err);
  sieve->reallocate();
#ifdef USE_DMCOMPLEX_ON
  /* DMPlex */
  err = DMPlexSymmetrize(newMesh);CHECK_PETSC_ERROR(err);
  err = DMPlexStratify(newMesh);CHECK_PETSC_ERROR(err);
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
  if (!faultSieveMesh->commRank()) {
    delete [] indices;
    err = PetscFree(indicesDM);CHECK_PETSC_ERROR(err);
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
  const ALE::Obj<SieveSubMesh::label_sequence>& fVertices2 = faultSieveMesh->depthStratum(0);
  assert(!fVertices2.isNull());
  SieveSubMesh::label_sequence::const_iterator fVertices2Begin = fVertices2->begin();
  SieveSubMesh::label_sequence::const_iterator fVertices2End   = fVertices2->end();

  PetscSection coordSection, newCoordSection;
  Vec          coordinatesVec, newCoordinatesVec;
  PetscScalar *coords, *newCoords;
  PetscInt     coordSize;
 
#ifdef USE_DMCOMPLEX_ON
  err = DMPlexGetCoordinateSection(complexMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(newMesh,     &newCoordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(complexMesh, &coordinatesVec);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(newCoordSection, vStart+extraCells, vEnd+extraCells+extraVertices);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetDof(newCoordSection, v+extraCells, dof);CHECK_PETSC_ERROR(err);
  }
#endif

  if (debug) coordinates->view("Coordinates without shadow vertices");
  for (SieveSubMesh::label_sequence::iterator v_iter = fVertices2Begin; v_iter != fVertices2End; ++v_iter) {
    coordinates->addPoint(vertexRenumber[*v_iter], coordinates->getFiberDimension(*v_iter));
    if (constraintCell) coordinates->addPoint(vertexLagrangeRenumber[*v_iter], coordinates->getFiberDimension(*v_iter));
#ifdef USE_DMCOMPLEX_ON
    PetscInt dof, v = *v_iter;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetDof(newCoordSection, vertexRenumberDM[v+faultVertexOffsetDM], dof);CHECK_PETSC_ERROR(err);
    if (constraintCell) {err = PetscSectionSetDof(newCoordSection, vertexLagrangeRenumberDM[v+faultVertexOffsetDM], dof);CHECK_PETSC_ERROR(err);}
#endif
  } // for
  sieveMesh->reallocate(coordinates);
#ifdef USE_DMCOMPLEX_ON
  err = PetscSectionSetUp(newCoordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(newCoordSection, &coordSize);CHECK_PETSC_ERROR(err);
  err = VecCreate(((PetscObject) newMesh)->comm, &newCoordinatesVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(newCoordinatesVec, coordSize, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(newCoordinatesVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordinatesVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecGetArray(newCoordinatesVec, &newCoords);CHECK_PETSC_ERROR(err);

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, newoff, d;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(newCoordSection, v+extraCells, &newoff);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      newCoords[newoff+d] = coords[off+d];
    }
  }
#endif
  SieveSubMesh::label_sequence::const_iterator fVertices2EndNew = fVertices2->end();
  for (SieveSubMesh::label_sequence::iterator v_iter = fVertices2Begin; v_iter != fVertices2EndNew; ++v_iter) {
    assert(coordinates->getFiberDimension(*v_iter) == coordinates->getFiberDimension(vertexRenumber[*v_iter]));
    coordinates->updatePoint(vertexRenumber[*v_iter], coordinates->restrictPoint(*v_iter));
    if (constraintCell) {
      assert(coordinates->getFiberDimension(*v_iter) == coordinates->getFiberDimension(vertexLagrangeRenumber[*v_iter]));
      coordinates->updatePoint(vertexLagrangeRenumber[*v_iter], coordinates->restrictPoint(*v_iter));
    } // if

#ifdef USE_DMCOMPLEX_ON
    PetscInt v = *v_iter, dof, off, newoff, d;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(newCoordSection, vertexRenumberDM[v+faultVertexOffsetDM], &newoff);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      newCoords[newoff+d] = coords[off+d];
    }
    if (constraintCell) {
      err = PetscSectionGetOffset(newCoordSection, vertexLagrangeRenumberDM[v+faultVertexOffsetDM], &newoff);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {
        newCoords[newoff+d] = coords[off+d];
      }
    } // if
#endif
  } // for
  err = VecRestoreArray(coordinatesVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(newCoordinatesVec, &newCoords);CHECK_PETSC_ERROR(err);
  err = DMSetCoordinatesLocal(newMesh, newCoordinatesVec);CHECK_PETSC_ERROR(err);
  err = VecDestroy(&newCoordinatesVec);CHECK_PETSC_ERROR(err);
  if (debug) coordinates->view("Coordinates with shadow vertices");

  logger.stagePop();

  mesh->setDMMesh(newMesh);
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

  err = DMPlexGetLabel(dm, labelName, &label);CHECK_PETSC_ERROR(err);
  // Completes the set of cells scheduled to be replaced
  //   Have to do internal fault vertices before fault boundary vertices, and this is the only thing I use faultBoundary for
  err = DMLabelCohesiveComplete(dm, label);CHECK_PETSC_ERROR(err);
  err = DMPlexConstructCohesiveCells(dm, label, &sdm);CHECK_PETSC_ERROR(err);
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
  assert(0 != faultMesh);

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultCreation");
  PetscErrorCode err;

  faultMesh->coordsys(mesh);

  DM dmMesh = mesh.dmMesh();
  assert(dmMesh);
#define OLD_PLEX
#ifndef OLD_PLEX
  DM              dmFaultMesh;
  {
    DMLabel         faultLabel;
    IS              cohesiveIS;
    const PetscInt *cohesiveCells;
    PetscInt        numCohesiveCells = 0, c;
    const char     *labelName = "_internal_fault";

    err = DMPlexGetStratumIS(dmMesh, "material-id", materialId, &cohesiveIS);CHECK_PETSC_ERROR(err);
    err = DMPlexCreateLabel(dmMesh, labelName);CHECK_PETSC_ERROR(err);
    err = DMPlexGetLabel(dmMesh, labelName, &faultLabel);CHECK_PETSC_ERROR(err);
    if (cohesiveIS) {
      err = ISGetSize(cohesiveIS, &numCohesiveCells);CHECK_PETSC_ERROR(err);
      err = ISGetIndices(cohesiveIS, &cohesiveCells);CHECK_PETSC_ERROR(err);
    }
    /* This is for uninterpolated meshes. For interpolated meshes, we just mark one face and its closure
       I am using the negative side vertices all the time, and not the Lagrange vertices. Does this matter???
    */
    for (c = 0; c < numCohesiveCells; ++c) {
      const PetscInt  point = cohesiveCells[c];
      const PetscInt *cone;
      PetscInt        coneSize, c, faceSize;

      err = DMPlexGetConeSize(dmMesh, point, &coneSize);CHECK_PETSC_ERROR(err);
      err = DMPlexGetCone(dmMesh, point, &cone);CHECK_PETSC_ERROR(err);
      if (!constraintCell) {
        faceSize = coneSize / 2;
      } else {
        faceSize = coneSize / 3;
      }
      for(c = 0; c < faceSize; ++c) {
        err = DMLabelSetValue(faultLabel, cone[c], 1);CHECK_PETSC_ERROR(err);
      }
    }
    if (cohesiveIS) {err = ISRestoreIndices(cohesiveIS, &cohesiveCells);CHECK_PETSC_ERROR(err);}
    err = ISDestroy(&cohesiveIS);CHECK_PETSC_ERROR(err);
    err = DMPlexCreateSubmesh(dmMesh, labelName, &dmFaultMesh);CHECK_PETSC_ERROR(err);
    err = DMPlexRemoveLabel(dmMesh, labelName, &faultLabel);CHECK_PETSC_ERROR(err);
    err = DMLabelDestroy(&faultLabel);CHECK_PETSC_ERROR(err);
    faultMesh->setDMMesh(dmFaultMesh);
  }
#endif

  const int debug = mesh.debug();
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  faultSieveMesh = new SieveSubMesh(mesh.comm(), mesh.dimension()-1, debug);
  const ALE::Obj<SieveMesh::sieve_type> ifaultSieve =
    new SieveMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!ifaultSieve.isNull());
  ALE::Obj<SieveFlexMesh> fault = 
    new SieveFlexMesh(mesh.comm(), mesh.dimension()-1, debug);
  assert(!fault.isNull());
  ALE::Obj<SieveFlexMesh::sieve_type> faultSieve =
    new SieveFlexMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!faultSieve.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& cohesiveCells =
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cohesiveCells.isNull());
  const SieveMesh::label_sequence::iterator cBegin = cohesiveCells->begin();
  const SieveMesh::label_sequence::iterator cEnd = cohesiveCells->end();
  const int sieveEnd = *std::max_element(sieve->getChart().begin(), sieve->getChart().end())+1;
  const int numFaces = cohesiveCells->size();
  int globalSieveEnd = 0;
  int globalFaceOffset = 0;

  // TODO: For multiple faults, this produces duplicate names. Not sure if we need to worry
  MPI_Allreduce((void *) &sieveEnd, (void *) &globalSieveEnd, 1,
		MPI_INT, MPI_SUM, sieve->comm());
  MPI_Scan((void *) &numFaces, (void *) &globalFaceOffset, 1,
	   MPI_INT, MPI_SUM, sieve->comm());
  int face = globalSieveEnd + globalFaceOffset - numFaces;

  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> cV(std::max(sieve->getMaxConeSize(), 1));

  for(SieveMesh::label_sequence::iterator c_iter = cBegin;
      c_iter != cEnd;
      ++c_iter) {
    sieve->cone(*c_iter, cV);
    const int coneSize = cV.getSize();
    const SieveMesh::point_type *cone = cV.getPoints();
    int color = 0;

    if (!constraintCell) {
      const int faceSize = coneSize / 2;
      assert(0 == coneSize % faceSize);

      // Use first vertices (negative side of the fault) for fault mesh
      for (int i = 0; i < faceSize; ++i)
        faultSieve->addArrow(cone[i], face, color++);
    } else {
      const int faceSize = coneSize / 3;
      assert(0 == coneSize % faceSize);

      // Use last vertices (contraints) for fault mesh
      for (int i = 2*faceSize; i < 3*faceSize; ++i)
        faultSieve->addArrow(cone[i], face, color++);
    } // if/else
    ++face;
    cV.clear();
  } // for
  fault->setSieve(faultSieve);
  logger.stagePop();
  logger.stagePush("FaultStratification");
  fault->stratify();
  logger.stagePop();
  logger.stagePush("FaultCreation");

  // Convert fault to an IMesh
  //   In general, renumbering[global point number] = local point number
  //   fRenumbering[mesh point] = fault mesh point
  SieveSubMesh::renumbering_type& fRenumbering =
    faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator fRenumberingEnd = 
    fRenumbering.end();
  faultSieveMesh->setSieve(ifaultSieve);
  //ALE::ISieveConverter::convertMesh(*fault, *faultSieveMesh, fRenumbering, true);
  {
    ALE::ISieveConverter::convertSieve(*fault->getSieve(),
				       *faultSieveMesh->getSieve(),
				       fRenumbering, true);
    logger.stagePop();
    logger.stagePush("FaultStratification");
    faultSieveMesh->stratify();
    logger.stagePop();
    logger.stagePush("FaultCreation");
    ALE::ISieveConverter::convertOrientation(*fault->getSieve(),
					     *faultSieveMesh->getSieve(),
					     fRenumbering,
					     fault->getArrowSection("orientation").ptr());
  }
#ifdef OLD_PLEX
  SieveMesh::renumbering_type convertRenumbering;
  DM dmFaultMesh;

  ALE::ISieveConverter::convertMesh(*fault, &dmFaultMesh, convertRenumbering, true);
  faultMesh->setDMMesh(dmFaultMesh);
#endif
  fault      = NULL;
  faultSieve = NULL;

  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  SieveSubMesh::label_sequence::iterator f_iter = faultCells->begin();

  // Update coordinates
  const ALE::Obj<topology::Mesh::RealSection>& coordinates =
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<topology::Mesh::RealSection>& fCoordinates =
    faultSieveMesh->getRealSection("coordinates");
  assert(!fCoordinates.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices =
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator vBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator vEnd = vertices->end();

  const ALE::Obj<SieveMesh::label_sequence>& fvertices =
    faultSieveMesh->depthStratum(0);
  assert(!fvertices.isNull());
  const point_type fvMin = *std::min_element(fvertices->begin(), 
					     fvertices->end());
  const point_type fvMax = *std::max_element(fvertices->begin(),
					     fvertices->end());
  
  fCoordinates->setChart(topology::Mesh::RealSection::chart_type(fvMin, 
								 fvMax+1));
  for (SieveMesh::label_sequence::iterator v_iter = vBegin;
       v_iter != vEnd;
       ++v_iter) {
    if (fRenumbering.find(*v_iter) == fRenumberingEnd)
      continue;
    fCoordinates->setFiberDimension(fRenumbering[*v_iter],
				    coordinates->getFiberDimension(*v_iter));
  } // for
  fCoordinates->allocatePoint();
  for(SieveMesh::label_sequence::iterator v_iter = vBegin;
      v_iter != vEnd;
      ++v_iter) {
    if (fRenumbering.find(*v_iter) == fRenumberingEnd)
      continue;
    fCoordinates->updatePoint(fRenumbering[*v_iter], 
			      coordinates->restrictPoint(*v_iter));
  }
#ifdef OLD_PLEX
  //faultSieveMesh->view("Parallel fault mesh");
  PetscInt numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices;
  mesh.getPointTypeSizes(&numNormalCells, &numCohesiveCells, &numNormalVertices, &numShadowVertices, &numLagrangeVertices);

  const SieveMesh::renumbering_type::const_iterator convertRenumberingEnd = convertRenumbering.end();
  PetscSection   coordSection, faultCoordSection;
  Vec            coordinateVec, faultCoordinateVec;
  PetscScalar   *a, *fa;
  PetscInt       pvStart, pvEnd, fvStart, fvEnd, n;

  err = DMPlexGetDepthStratum(dmMesh, 0, &pvStart, &pvEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmFaultMesh, 0, &fvStart, &fvEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMPlexGetCoordinateSection(dmFaultMesh, &faultCoordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionSetChart(faultCoordSection, fvStart, fvEnd);CHECK_PETSC_ERROR(err);
  for(PetscInt v = pvStart; v < pvEnd; ++v) {
    PetscInt sievePoint = convertDMToSievePointNumbering(v, numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices);
    PetscInt dof;

    if (convertRenumbering.find(sievePoint) == convertRenumberingEnd) continue;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionSetDof(faultCoordSection, convertRenumbering[sievePoint], dof);CHECK_PETSC_ERROR(err);
  }
  err = PetscSectionSetUp(faultCoordSection);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetStorageSize(faultCoordSection, &n);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordinateVec);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmFaultMesh, &faultCoordinateVec);CHECK_PETSC_ERROR(err);
  err = VecSetSizes(faultCoordinateVec, n, PETSC_DETERMINE);CHECK_PETSC_ERROR(err);
  err = VecSetFromOptions(faultCoordinateVec);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordinateVec, &a);CHECK_PETSC_ERROR(err);
  err = VecGetArray(faultCoordinateVec, &fa);CHECK_PETSC_ERROR(err);
  for(PetscInt v = pvStart; v < pvEnd; ++v) {
    PetscInt sievePoint = convertDMToSievePointNumbering(v, numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices);
    PetscInt dof, off, foff;

    if (convertRenumbering.find(sievePoint) == convertRenumberingEnd) continue;
    err = PetscSectionGetDof(coordSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(faultCoordSection, convertRenumbering[sievePoint], &foff);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {
      fa[foff+d] = a[off+d];
    }
  }
  err = VecRestoreArray(coordinateVec, &a);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(faultCoordinateVec, &fa);CHECK_PETSC_ERROR(err);
#endif

#ifdef OLD_PLEX
  // Have to make subpointMap here: renumbering[original] = fault
  DMLabel   subpointMap;
  PetscInt *renum;
  PetscInt  pStart, pEnd, fcStart, fcEnd;

  err = DMLabelCreate("subpoint_map", &subpointMap);CHECK_PETSC_ERROR(err);
  err = DMPlexGetChart(dmFaultMesh, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetDepthStratum(dmFaultMesh, 0, &fvStart, &fvEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexGetHeightStratum(dmFaultMesh, 0, &fcStart, &fcEnd);CHECK_PETSC_ERROR(err);
  assert(convertRenumbering.size() == pEnd-pStart);
  err = PetscMalloc((pEnd-pStart) * sizeof(PetscInt), &renum);CHECK_PETSC_ERROR(err);
#if 0
  mesh.sieveMesh()->view("Sieve Mesh");
  err = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL);CHECK_PETSC_ERROR(err);
  err = DMView(mesh.dmMesh(), PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
  err = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHECK_PETSC_ERROR(err);
#endif
  for(SieveMesh::renumbering_type::const_iterator p_iter = convertRenumbering.begin(); p_iter != convertRenumbering.end(); ++p_iter) {
    const PetscInt sievePoint = p_iter->first;
    const PetscInt faultPoint = p_iter->second;
    PetscInt       dmPoint    = convertSieveToDMPointNumbering(sievePoint, numNormalCells, numCohesiveCells, numNormalVertices, numShadowVertices, numLagrangeVertices);
    renum[faultPoint] = dmPoint;
    //std::cout << "renum["<<faultPoint<<"]: "<<dmPoint<<std::endl;
  }
  for(PetscInt p = 1; p < pEnd-pStart; ++p) {
    if (renum[p-1] == -1) continue;
    assert(renum[p] > renum[p-1]);
  }
  for(PetscInt p = fcStart; p < fcEnd; ++p) {
    err = DMLabelSetValue(subpointMap, renum[p], mesh.dimension());CHECK_PETSC_ERROR(err);
  }
  for(PetscInt p = fvStart; p < fvEnd; ++p) {
    err = DMLabelSetValue(subpointMap, renum[p], 0);CHECK_PETSC_ERROR(err);
  }
  err = PetscFree(renum);CHECK_PETSC_ERROR(err);
  err = DMPlexSetSubpointMap(dmFaultMesh, subpointMap);CHECK_PETSC_ERROR(err);
#endif

  // Update dimensioned coordinates if they exist.
  if (sieveMesh->hasRealSection("coordinates_dimensioned")) {
    const ALE::Obj<topology::Mesh::RealSection>& coordinatesDim =
      sieveMesh->getRealSection("coordinates_dimensioned");
    assert(!coordinatesDim.isNull());
    const ALE::Obj<topology::Mesh::RealSection>& fCoordinatesDim =
      faultSieveMesh->getRealSection("coordinates_dimensioned");
    assert(!fCoordinatesDim.isNull());

    fCoordinatesDim->setChart(topology::Mesh::RealSection::chart_type(fvMin,
								      fvMax+1));
    for (SieveMesh::label_sequence::iterator v_iter = vBegin;
	 v_iter != vEnd;
	 ++v_iter) {
      if (fRenumbering.find(*v_iter) == fRenumberingEnd)
	continue;
      fCoordinatesDim->setFiberDimension(fRenumbering[*v_iter],
					 coordinatesDim->getFiberDimension(*v_iter));
    } // for
    fCoordinatesDim->allocatePoint();
    for(SieveMesh::label_sequence::iterator v_iter = vBegin;
	v_iter != vEnd;
	++v_iter) {
      if (fRenumbering.find(*v_iter) == fRenumberingEnd)
	continue;
      assert(fCoordinatesDim->getFiberDimension(fRenumbering[*v_iter]) ==
	     coordinatesDim->getFiberDimension(*v_iter));
      fCoordinatesDim->updatePoint(fRenumbering[*v_iter], 
				   coordinatesDim->restrictPoint(*v_iter));
    } // for
    //faultSieveMesh->view("Parallel fault mesh");
  } // if

  // Create the parallel overlap
  //   Can I figure this out in a nicer way?
  ALE::Obj<SieveSubMesh::send_overlap_type> sendParallelMeshOverlap =
    faultSieveMesh->getSendOverlap();
  assert(!sendParallelMeshOverlap.isNull());
  ALE::Obj<SieveSubMesh::recv_overlap_type> recvParallelMeshOverlap =
    faultSieveMesh->getRecvOverlap();
  assert(!recvParallelMeshOverlap.isNull());

  // Must process the renumbering local --> fault to global --> fault
  SieveMesh::renumbering_type& renumbering = sieveMesh->getRenumbering();
  SieveMesh::renumbering_type gRenumbering;

  if (renumbering.size()) {
    //std::cout << "Using renumbering to construct Fault Overlap" << std::endl;
    const SieveMesh::renumbering_type::const_iterator renumberingEnd =
      renumbering.end();
    for (SieveMesh::renumbering_type::const_iterator r_iter = renumbering.begin();
         r_iter != renumberingEnd;
         ++r_iter)
      if (fRenumbering.find(r_iter->second) != fRenumbering.end())
        gRenumbering[r_iter->first] = fRenumbering[r_iter->second];
  } else {
    //std::cout << "Using new numbering to construct Fault Overlap" << std::endl;
    const SieveMesh::sieve_type::chart_type& chart = sieveMesh->getSieve()->getChart();
    const ALE::Obj<SieveMesh::numbering_type>& globalNumbering = 
      sieveMesh->getFactory()->getNumbering(sieveMesh, -1);
    assert(!globalNumbering.isNull());
    for(SieveMesh::point_type p = chart.min(); p < chart.max(); ++p) {
      if (fRenumbering.find(p) != fRenumbering.end()) {
        gRenumbering[globalNumbering->getIndex(p)] = fRenumbering[p];
      } // if
    } // for
  } // if/else

  ALE::SetFromMap<SieveMesh::renumbering_type> globalPoints(gRenumbering);
  ALE::OverlapBuilder<>::constructOverlap(globalPoints, gRenumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  faultSieveMesh->setCalculatedOverlap(true);

#if 0 // Seems to break unit tests.
  // Consistency check for parallel overlap.
  if (fRenumbering.size() > 0) {
    if (renumbering.size() <= 0 ||
	gRenumbering.size() <= 0) {
      throw std::logic_error("Inconsistent data when computing overlap for "
			     "parallel fault mesh.");
    } // if
  } // if
#endif
  
#if 0 // DEBUGGING
  sendParallelMeshOverlap->view("Send parallel fault overlap");
  recvParallelMeshOverlap->view("Recv parallel fault overlap");
#endif

  logger.stagePop();
} // createFaultParallel


// End of file
