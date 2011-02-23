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
// Copyright (c) 2010 University of California, Davis
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
                                              ALE::Obj<FlexMesh>& faultBoundary,
                                              const topology::Mesh& mesh,
                                              const ALE::Obj<topology::Mesh::IntSection>& groupField,
					      const bool flipFault)
{ // createFault
  assert(0 != faultMesh);
  assert(!groupField.isNull());

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultCreation");

  faultMesh->coordsys(mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();
  faultSieveMesh =
    new SieveSubMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<SieveSubMesh::sieve_type> ifaultSieve = 
    new SieveMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!ifaultSieve.isNull());
  ALE::Obj<FlexMesh> fault =
    new FlexMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());
  assert(!fault.isNull());
  ALE::Obj<FlexMesh::sieve_type> faultSieve  =
    new FlexMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!faultSieve.isNull());
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
  logger.stagePush("FaultStratification");
  fault->stratify();
  logger.stagePop();
  logger.stagePush("FaultCreation");
  if (debug)
    fault->view("Fault mesh");

  faultBoundary = ALE::Selection<FlexMesh>::boundary(fault);
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

  logger.stagePop();
} // createFault

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::create(topology::Mesh* mesh,
                                         const topology::SubMesh& faultMesh,
                                         const ALE::Obj<FlexMesh>& faultBoundary,
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

  typedef ALE::SieveAlg<FlexMesh> sieveAlg;
  typedef ALE::Selection<FlexMesh> selection;

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("FaultCreation");

  const ALE::Obj<SieveMesh>& sieveMesh = mesh->sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  assert(!faultSieveMesh.isNull());  

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<SieveSubMesh::sieve_type> ifaultSieve = 
    faultSieveMesh->getSieve();
  assert(!ifaultSieve.isNull());

  const int  depth = sieveMesh->depth();
  assert(!sieveMesh->heightStratum(0).isNull());
  const int numCells = sieveMesh->heightStratum(0)->size();
  int numCorners = 0; // The number of vertices in a mesh cell
  int faceSize = 0; // The number of vertices in a mesh face
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int* indices = 0; // The indices of a face vertex set in a cell
  const int debug = mesh->debug();
  int oppositeVertex = 0;    // For simplices, the vertex opposite a given face
  TopologyOps::PointArray origVertices;
  TopologyOps::PointArray faceVertices;
  TopologyOps::PointArray neighborVertices;

  if (!faultSieveMesh->commRank()) {
    assert(!faultSieveMesh->heightStratum(1).isNull());
    const SieveSubMesh::point_type p = *faultSieveMesh->heightStratum(1)->begin();

    numCorners = sieveMesh->getNumCellCorners();
    faceSize = selection::numFaceVertices(sieveMesh);
    indices = new int[faceSize];
    numFaultCorners = faultSieveMesh->getNumCellCorners(p, faultSieveMesh->depth(p));
  }
  //faultSieveMesh->view("Serial fault mesh");

  // Add new shadow vertices and possibly Lagrange multipler vertices
  const ALE::Obj<SieveSubMesh::label_sequence>& fVertices       = faultSieveMesh->depthStratum(0);
  assert(!fVertices.isNull());
  const SieveSubMesh::label_sequence::const_iterator fVerticesBegin = 
    fVertices->begin();
  const SieveSubMesh::label_sequence::const_iterator fVerticesEnd = 
    fVertices->end();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const ALE::Obj<std::set<std::string> >& groupNames = 
    sieveMesh->getIntSections();
  assert(!groupNames.isNull());
  const int numFaultVertices = fVertices->size();
  std::map<point_type,point_type> vertexRenumber;
  std::map<point_type,point_type> vertexLagrangeRenumber;
  std::map<point_type,point_type> cellRenumber;
  if (firstFaultVertex == 0) {
    firstFaultVertex    += sieve->getBaseSize() + sieve->getCapSize();
    firstLagrangeVertex += firstFaultVertex;
    firstFaultCell      += firstFaultVertex;
  }

  for(SieveSubMesh::label_sequence::iterator v_iter = fVerticesBegin;
      v_iter != fVerticesEnd;
      ++v_iter, ++firstFaultVertex) {
    vertexRenumber[*v_iter] = firstFaultVertex;
    if (debug) 
      std::cout << "Duplicating " << *v_iter << " to "
		<< vertexRenumber[*v_iter] << std::endl;

    logger.stagePop();
    logger.stagePush("FaultStratification");
    // Add shadow and constraint vertices (if they exist) to group
    // associated with fault
    groupField->addPoint(firstFaultVertex, 1);
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultVertex, 1);
    sieveMesh->setDepth(firstFaultVertex, 0);
#endif
    if (constraintCell) {
      vertexLagrangeRenumber[*v_iter] = firstLagrangeVertex;
      groupField->addPoint(firstLagrangeVertex, 1);
#if defined(FAST_STRATIFY)
      // OPTIMIZATION
      sieveMesh->setHeight(firstLagrangeVertex, 1);
      sieveMesh->setDepth(firstLagrangeVertex, 0);
#endif
      ++firstLagrangeVertex;
    } // if
    logger.stagePop();
    logger.stagePush("FaultCreation");

    // Add shadow vertices to other groups, don't add constraint
    // vertices (if they exist) because we don't want BC, etc to act
    // on constraint vertices
    const std::set<std::string>::const_iterator namesEnd = groupNames->end();
    for(std::set<std::string>::const_iterator name = groupNames->begin();
       name != namesEnd;
	++name) {
      const ALE::Obj<IntSection>& group = sieveMesh->getIntSection(*name);
      assert(!group.isNull());
      if (group->getFiberDimension(*v_iter))
        group->addPoint(firstFaultVertex, 1);
    } // for
  } // for
  const std::set<std::string>::const_iterator namesEnd = groupNames->end();
  for(std::set<std::string>::const_iterator name = groupNames->begin();
      name != namesEnd;
      ++name) {
    sieveMesh->reallocate(sieveMesh->getIntSection(*name));
  } // for

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<SieveSubMesh::label_sequence>& faces =
    faultSieveMesh->heightStratum(1);
  assert(!faces.isNull());
  const SieveSubMesh::label_sequence::const_iterator facesBegin = faces->begin();
  const SieveSubMesh::label_sequence::const_iterator facesEnd = faces->end();
  const ALE::Obj<SieveFlexMesh::label_type>& material = 
    sieveMesh->getLabel("material-id");
  assert(!material.isNull());
  const int firstCohesiveCell = firstFaultCell;
  TopologyOps::PointSet replaceCells;
  TopologyOps::PointSet noReplaceCells;
  TopologyOps::PointSet replaceVertices;
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV2(std::max(1, ifaultSieve->getMaxSupportSize()));
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> cV2(*ifaultSieve, (size_t) pow(std::max(1, ifaultSieve->getMaxConeSize()), faultSieveMesh->depth()));
  std::set<SieveFlexMesh::point_type> faceSet;

  for(SieveSubMesh::label_sequence::iterator f_iter = facesBegin;
      f_iter != facesEnd;
      ++f_iter, ++firstFaultCell) {
    const point_type face = *f_iter;
    if (debug)
      std::cout << "Considering fault face " << face << std::endl;
    ifaultSieve->support(face, sV2);
    const point_type *cells = sV2.getPoints();
    point_type cell = cells[0];
    point_type otherCell = cells[1];

    if (debug)
      std::cout << "  Checking orientation against cell " << cell << std::endl;
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*ifaultSieve,
								 face, cV2);
    const int coneSize = cV2.getSize();
    const point_type *faceCone = cV2.getPoints();
    //ifaultSieve->cone(face, cV2);
    //const int coneSize = cV2.getSize() ? cV2.getSize() : 1;
    //const point_type *faceCone = cV2.getSize() ? cV2.getPoints() : &face;
    bool found = true;

    for(int i = 0; i < coneSize; ++i)
      faceSet.insert(faceCone[i]);
    selection::getOrientedFace(sieveMesh, cell, &faceSet, numCorners, indices,
			       &origVertices, &faceVertices);
    if (faceVertices.size() != coneSize) {
      std::cout << "Invalid size for faceVertices " << faceVertices.size()
		<< " of face " << face << "should be " << coneSize << std::endl;
      std::cout << "  firstCohesiveCell " << firstCohesiveCell << " firstFaultCell " 
		<< firstFaultCell << " numFaces " << faces->size() << std::endl;
      std::cout << "  faceSet:" << std::endl;
      for(std::set<SieveFlexMesh::point_type>::const_iterator p_iter = faceSet.begin();
	  p_iter != faceSet.end();
	  ++p_iter) {
        std::cout << "    " << *p_iter << std::endl;
      } // if
      std::cout << "  cell cone:" << std::endl;
      ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> cV(std::max(1, sieve->getMaxConeSize()));
      sieve->cone(cell, cV);
      const int coneSize2 = cV.getSize();
      const point_type *cellCone  = cV.getPoints();

      for(int c = 0; c < coneSize2; ++c)
        std::cout << "    " << cellCone[c] << std::endl;
      std::cout << "  fault cell support:" << std::endl;
      ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV(std::max(1, ifaultSieve->getMaxSupportSize()));
      ifaultSieve->support(face, sV);
      const int supportSize2 = sV.getSize();
      const point_type *cellSupport  = sV.getPoints();
      for(int s = 0; s < supportSize2; ++s)
        std::cout << "    " << cellSupport[s] << std::endl;
    } // if
    assert(faceVertices.size() == coneSize);
    faceSet.clear();
    //selection::getOrientedFace(sieveMesh, cell, &vertexRenumber, numCorners,
    //			       indices, &origVertices, &faceVertices);

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
        if (debug) std::cout << "    Checking " << faceCone[c] << " against "
			     << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != faceCone[c]) {
          found = false;
          break;
        } // if
      } // for
    } // if/else

    if (found) {
      if (debug)
	std::cout << "  Choosing other cell" << std::endl;
      point_type tmpCell = otherCell;
      otherCell = cell;
      cell = tmpCell;
    } else {
      if (debug)
	std::cout << "  Verifing reverse orientation" << std::endl;
      found = true;
      int v = 0;
      if (numFaultCorners > 0) {
        // Locate first vertex
        while((v < numFaultCorners) && (faceVertices[v] != faceCone[coneSize-1]))
	  ++v;
        for(int c = coneSize-1; c >= 0; --c, ++v) {
          if (debug)
	    std::cout << "    Checking " << faceCone[c] << " against "
		      << faceVertices[v%numFaultCorners] << std::endl;
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
    cellRenumber[cell] = firstFaultCell;
    // Adding cohesive cell (not interpolated)
    if (debug)
      std::cout << "  Creating cohesive cell " << firstFaultCell << std::endl;
    for (int c = 0; c < coneSize; ++c) {
      if (debug)
	std::cout << "    vertex " << faceCone[c] << std::endl;
      sieve->addArrow(faceCone[c], firstFaultCell);
    } // for
    for (int c = 0; c < coneSize; ++c) {
      if (debug)
	std::cout << "    shadow vertex " << vertexRenumber[faceCone[c]] << std::endl;
      sieve->addArrow(vertexRenumber[faceCone[c]], firstFaultCell, true);
    } // for
    if (constraintCell) {
      for (int c = 0; c < coneSize; ++c) {
        if (debug)
	  std::cout << "    Lagrange vertex " << vertexLagrangeRenumber[faceCone[c]] << std::endl;
        sieve->addArrow(vertexLagrangeRenumber[faceCone[c]], firstFaultCell, true);
      } // for
    } // if
    // TODO: Need to reform the material label when sieve is reallocated
    sieveMesh->setValue(material, firstFaultCell, materialId);
    logger.stagePop();
    logger.stagePush("FaultStratification");
#if defined(FAST_STRATIFY)
    // OPTIMIZATION
    sieveMesh->setHeight(firstFaultCell, 0);
    sieveMesh->setDepth(firstFaultCell, 1);
#endif
    logger.stagePop();
    logger.stagePush("FaultCreation");
    sV2.clear();
    cV2.clear();
  } // for
  // This completes the set of cells scheduled to be replaced
  TopologyOps::PointSet replaceCellsBase(replaceCells);

  const ALE::Obj<FlexMesh::label_sequence>& faultBdVerts =
    faultBoundary->depthStratum(0);
  assert(!faultBdVerts.isNull());
  TopologyOps::PointSet faultBdVertices;

  faultBdVertices.insert(faultBdVerts->begin(), faultBdVerts->end());
  TopologyOps::PointSet::const_iterator rVerticesEnd = replaceVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin();
      v_iter != rVerticesEnd; ++v_iter) {
    if (faultBdVertices.find(*v_iter) != faultBdVertices.end())
      continue;
    TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize,
			       firstCohesiveCell, replaceCells, noReplaceCells,
			       debug);
  } // for
  const TopologyOps::PointSet::const_iterator fbdVerticesEnd = 
    faultBdVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter=faultBdVertices.begin();
      v_iter != fbdVerticesEnd;
      ++v_iter) {
    TopologyOps::classifyCells(sieve, *v_iter, depth, faceSize,
			       firstCohesiveCell, replaceCells, noReplaceCells,
			       debug);
  } // for
  // Add new arrows for support of replaced vertices
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV(std::max(1, sieve->getMaxSupportSize()));

  rVerticesEnd = replaceVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin();
      v_iter != rVerticesEnd;
      ++v_iter) {
    sieve->support(*v_iter, sV);
    const point_type *support = sV.getPoints();

    if (debug)
      std::cout << "  Checking support of " << *v_iter << std::endl;
    const int sVSize = sV.getSize();
    for (int s = 0; s < sVSize; ++s) {
      if (replaceCells.find(support[s]) != replaceCells.end()) {
        if (debug)
	  std::cout << "    Adding new support " << vertexRenumber[*v_iter]
		    << " --> " << support[s] << std::endl;
        sieve->addArrow(vertexRenumber[*v_iter], support[s], true);
      } else {
        if (debug)
	  std::cout << "    Keeping same support " << *v_iter<<","
		    << vertexRenumber[*v_iter] << " --> " << support[s]
		    << std::endl;
      } // if/else
    } // for
    sV.clear();
  }
  sieve->reallocate();

  // More checking
  const bool firstFault = !sieveMesh->hasRealSection("replaced_cells");
  const ALE::Obj<topology::Mesh::RealSection>& replacedCells = 
    sieveMesh->getRealSection("replaced_cells");
  assert(!replacedCells.isNull());
  TopologyOps::PointSet cellNeighbors;
	 
  if (firstFault) {
    const ALE::Obj<SieveMesh::label_sequence>& cells = 
      sieveMesh->heightStratum(0);
    assert(!cells.isNull());

    replacedCells->setChart(topology::Mesh::RealSection::chart_type(*std::min_element(cells->begin(), cells->end()), *std::max_element(cells->begin(), cells->end())+1));
    replacedCells->setFiberDimension(cells, 1);
    replacedCells->allocatePoint();
  } // if
	 
  const TopologyOps::PointSet::const_iterator noRCellsEnd = noReplaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = noReplaceCells.begin();
      c_iter != noRCellsEnd;
      ++c_iter) {
    const double minusOne = -1.0;
    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &minusOne);
    } else {
      const double minusTwo = -2.0;
      replacedCells->updatePoint(*c_iter, &minusTwo);
    } // if/else
  } // for

  TopologyOps::PointSet::const_iterator rCellsEnd = replaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = replaceCells.begin();
      c_iter != rCellsEnd;
      ++c_iter) {
    if (replaceCellsBase.find(*c_iter) != replaceCellsBase.end()) {
      const double one = 1.0;
      if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
        replacedCells->updatePoint(*c_iter, &one);
      } else {
        const double two = 2.0;
        replacedCells->updatePoint(*c_iter, &two);
      } // if/else
      continue;
    } // if
    const double ten = 10.0;
    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &ten);
    } else {
      const double twenty = 20.0;
      replacedCells->updatePoint(*c_iter, &twenty);
    } // if/else
    // There should be a way to check for boundary elements
    if (mesh->dimension() == 1) {
      if (cellNeighbors.size() > 2) {
        std::cout << "Cell " << *c_iter
		  << " has an invalid number of neighbors "
		  << cellNeighbors.size() << std::endl;
        throw ALE::Exception("Invalid number of neighbors");
      } // if
    } else if (mesh->dimension() == 2) {
      if (numCorners == 3) {
        if (cellNeighbors.size() > 3) {
          std::cout << "Cell " << *c_iter
		    << " has an invalid number of neighbors "
		    << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
	} // if
      } else if (numCorners == 4 || numCorners == 9) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter
		    << " has an invalid number of neighbors "
		    << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } // if/else
    } else if (mesh->dimension() == 3) {
      if (numCorners == 4) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter
		    << " has an invalid number of neighbors "
		    << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } else if (numCorners == 8) {
        if (cellNeighbors.size() > 6) {
          std::cout << "Cell " << *c_iter
		    << " has an invalid number of neighbors "
		    << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        } // if
      } // if/else
    } // if/else
  } // for
  ReplaceVisitor<SieveMesh::sieve_type,std::map<SieveMesh::point_type,SieveMesh::point_type> > rVc(vertexRenumber, std::max(1, sieve->getMaxConeSize()), debug);
  
  rCellsEnd = replaceCells.end();
  for (TopologyOps::PointSet::const_iterator c_iter = replaceCells.begin();
       c_iter != rCellsEnd;
       ++c_iter) {
    sieve->cone(*c_iter, rVc);
    if (rVc.mappedPoint()) {
      if (debug)
	std::cout << "  Replacing cell " << *c_iter << std::endl;
      sieve->setCone(rVc.getPoints(), *c_iter);
    } // if
    rVc.clear();
  } // for
  ReplaceVisitor<SieveMesh::sieve_type,std::map<SieveMesh::point_type,SieveMesh::point_type> > rVs(cellRenumber, std::max(1, sieve->getMaxSupportSize()), debug);

  rVerticesEnd = replaceVertices.end();
  for (TopologyOps::PointSet::const_iterator v_iter = replaceVertices.begin();
       v_iter != rVerticesEnd;
       ++v_iter) {
    sieve->support(*v_iter, rVs);
    if (rVs.mappedPoint()) {
      if (debug)
	std::cout << "  Replacing support for " << *v_iter << std::endl;
      sieve->setSupport(*v_iter, rVs.getPoints());
    } else {
      if (debug)
	std::cout << "  Not replacing support for " << *v_iter << std::endl;
    } // if/else
    rVs.clear();
  } // for
  if (!faultSieveMesh->commRank())
    delete [] indices;
#if !defined(FAST_STRATIFY)
  logger.stagePop();
  logger.stagePush("FaultStratification");
  sieveMesh->stratify();
  logger.stagePop();
  logger.stagePush("FaultCreation");
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
    for (std::map<int,int>::const_iterator v_iter = vertexRenumber.begin();
	 v_iter != vRenumberEnd;
	 ++v_iter)
      sieveMesh->setValue(label, v_iter->second, 0);
  } // if/else
  if (debug)
    mesh->view("Mesh with Cohesive Elements");

  // Fix coordinates
  const ALE::Obj<topology::Mesh::RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& fVertices2 =
    faultSieveMesh->depthStratum(0);
  assert(!fVertices2.isNull());
  SieveSubMesh::label_sequence::const_iterator fVertices2Begin = 
    fVertices2->begin();
  SieveSubMesh::label_sequence::const_iterator fVertices2End = 
    fVertices2->end();

  if (debug)
    coordinates->view("Coordinates without shadow vertices");
  for (SieveSubMesh::label_sequence::iterator v_iter = fVertices2Begin;
      v_iter != fVertices2End;
      ++v_iter) {
    coordinates->addPoint(vertexRenumber[*v_iter],
			  coordinates->getFiberDimension(*v_iter));
    if (constraintCell)
      coordinates->addPoint(vertexLagrangeRenumber[*v_iter],
			    coordinates->getFiberDimension(*v_iter));
  } // for
  sieveMesh->reallocate(coordinates);
  SieveSubMesh::label_sequence::const_iterator fVertices2EndNew = 
    fVertices2->end();
  for (SieveSubMesh::label_sequence::iterator v_iter = fVertices2Begin;
       v_iter != fVertices2EndNew;
       ++v_iter) {
    coordinates->updatePoint(vertexRenumber[*v_iter], 
			     coordinates->restrictPoint(*v_iter));
    if (constraintCell)
      coordinates->updatePoint(vertexLagrangeRenumber[*v_iter],
			       coordinates->restrictPoint(*v_iter));
  } // for
  if (debug)
    coordinates->view("Coordinates with shadow vertices");

  logger.stagePop();
} // create

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

  faultMesh->coordsys(mesh);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh->sieveMesh();

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  faultSieveMesh = 
    new SieveSubMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());
  const ALE::Obj<SieveMesh::sieve_type> ifaultSieve =
    new SieveMesh::sieve_type(sieve->comm(), sieve->debug());
  assert(!ifaultSieve.isNull());
  ALE::Obj<FlexMesh> fault = 
    new FlexMesh(mesh.comm(), mesh.dimension()-1, mesh.debug());
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
  //faultSieveMesh->view("Parallel fault mesh");

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

  const SieveMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();
  for (SieveMesh::renumbering_type::const_iterator r_iter = renumbering.begin();
       r_iter != renumberingEnd;
       ++r_iter)
    if (fRenumbering.find(r_iter->second) != fRenumbering.end())
      gRenumbering[r_iter->first] = fRenumbering[r_iter->second];

  ALE::SetFromMap<SieveMesh::renumbering_type> globalPoints(gRenumbering);
  ALE::OverlapBuilder<>::constructOverlap(globalPoints, gRenumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  faultSieveMesh->setCalculatedOverlap(true);
  //sendParallelMeshOverlap->view("Send parallel fault overlap");
  //recvParallelMeshOverlap->view("Recv parallel fault overlap");

  logger.stagePop();
} // createFaultParallel


// End of file
