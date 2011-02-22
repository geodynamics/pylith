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

#include "TopologyOps.hh" // implementation of object methods

#include "TopologyVisitors.hh" // USES ClassifyVisitor

#include <Selection.hh> // Algorithms for submeshes

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
template<class InputPoints>
bool
pylith::faults::TopologyOps::compatibleOrientation(const ALE::Obj<SieveMesh>& mesh,
							const SieveMesh::point_type& p,
							const SieveMesh::point_type& q,
							const int numFaultCorners,
							const int faultFaceSize,
							const int faultDepth,
						   const ALE::Obj<InputPoints>& points,
							int indices[],
							PointArray *origVertices,
							PointArray *faceVertices,
							PointArray *neighborVertices)
{
  typedef ALE::Selection<SieveMesh> selection;
  const int debug = mesh->debug();
  bool compatible;

  bool eOrient = selection::getOrientedFace(mesh, p, points, numFaultCorners, indices, origVertices, faceVertices);
  bool nOrient = selection::getOrientedFace(mesh, q, points, numFaultCorners, indices, origVertices, neighborVertices);

  if (faultFaceSize > 1) {
    if (debug) {
      for(PointArray::iterator v_iter = faceVertices->begin(); v_iter != faceVertices->end(); ++v_iter) {
        std::cout << "  face vertex " << *v_iter << std::endl;
      }
      for(PointArray::iterator v_iter = neighborVertices->begin(); v_iter != neighborVertices->end(); ++v_iter) {
        std::cout << "  neighbor vertex " << *v_iter << std::endl;
      }
    }
    compatible = !(*faceVertices->begin() == *neighborVertices->begin());
  } else {
    compatible = !(nOrient == eOrient);
  }
  return compatible;
}

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::computeCensoredDepth(const ALE::Obj<SieveMesh::label_type>& depth,
						       const ALE::Obj<SieveMesh::sieve_type>& sieve,
						       const SieveMesh::point_type& firstCohesiveCell)
{
  SieveMesh::DepthVisitor d(*sieve, firstCohesiveCell, *depth);

  sieve->roots(d);
  while(d.isModified()) {
    // FIX: Avoid the copy here somehow by fixing the traversal
    std::vector<SieveMesh::point_type> modifiedPoints(d.getModifiedPoints().begin(), d.getModifiedPoints().end());

    d.clear();
    sieve->support(modifiedPoints, d);
  }
}

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::classifyCells(const ALE::Obj<SieveMesh::sieve_type>& sieve,
                                                const SieveMesh::point_type& vertex,
                                                const int depth,
                                                const int faceSize,
                                                const SieveMesh::point_type& firstCohesiveCell,
                                                PointSet& replaceCells,
                                                PointSet& noReplaceCells,
                                                const int debug)
{
  // Replace all cells on a given side of the fault with a vertex on the fault
  ClassifyVisitor<SieveMesh::sieve_type> cV(*sieve, replaceCells, noReplaceCells,
					    firstCohesiveCell, faceSize, debug);
  const PointSet& vReplaceCells   = cV.getReplaceCells();
  const PointSet& vNoReplaceCells = cV.getNoReplaceCells();

  if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
  sieve->support(vertex, cV);
  cV.setMode(false);
  const int classifyTotal = cV.getSize();
  int       classifySize  = vReplaceCells.size() + vNoReplaceCells.size();

  while(cV.getModified() && (classifySize < classifyTotal)) {
    cV.reset();
    sieve->support(vertex, cV);
    if (debug) {
      std::cout << "classifySize: " << classifySize << std::endl;
      std::cout << "classifyTotal: " << classifyTotal << std::endl;
      std::cout << "vReplaceCells.size: " << vReplaceCells.size() << std::endl;
      std::cout << "vNoReplaceCells.size: " << vNoReplaceCells.size() << std::endl;
    }
    assert(classifySize < vReplaceCells.size() + vNoReplaceCells.size());
    classifySize = vReplaceCells.size() + vNoReplaceCells.size();
    assert(classifySize <= classifyTotal);
  }
  replaceCells.insert(vReplaceCells.begin(), vReplaceCells.end());
  // More checking
  noReplaceCells.insert(vNoReplaceCells.begin(), vNoReplaceCells.end());
}

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::createFaultSieveFromVertices(const int dim,
                                                               const int firstCell,
                                                               const PointSet& faultVertices,
                                                               const ALE::Obj<SieveMesh>& mesh,
                                                               const ALE::Obj<FlexMesh::arrow_section_type>& orientation,
                                                               const ALE::Obj<FlexMesh::sieve_type>& faultSieve,
							       const bool flipFault)
{
  typedef ALE::Selection<FlexMesh> selection;
  const ALE::Obj<SieveMesh::sieve_type>& sieve = mesh->getSieve();
  const PointSet::const_iterator fvBegin    = faultVertices.begin();
  const PointSet::const_iterator fvEnd      = faultVertices.end();
  int                            curCell    = firstCell;
  int                            curVertex  = 0;
  int                            newElement = curCell + dim*faultVertices.size();
  int                            o          = 1;
  FlexMesh::point_type          f          = firstCell;
  const int                      debug      = mesh->debug();
  ALE::Obj<PointSet>                  face       = new PointSet();
  int                            numCorners = 0;    // The number of vertices in a mesh cell
  int                            faceSize   = 0;    // The number of vertices in a mesh face
  int                           *indices    = NULL; // The indices of a face vertex set in a cell
  std::map<int,int*>             curElement;
  std::map<int,PointArray>       bdVertices;
  std::map<int,PointArray>       faultFaces;
  std::map<int,oPointArray>      oFaultFaces;
  PointSet                       faultCells;
  PointArray                     origVertices;
  PointArray                     faceVertices;

  //faultSieve->setDebug(2);
  if (!faultSieve->commRank()) {
    numCorners = mesh->getNumCellCorners();
    faceSize   = selection::numFaceVertices(mesh);
    indices    = new int[faceSize];
  }

  curElement[0]   = &curVertex;
  curElement[dim] = &curCell;
  for(int d = 1; d < dim; d++) {
    curElement[d] = &newElement;
  }

  // This only works for uninterpolated meshes
  assert((mesh->depth() == 1) || (mesh->depth() == -1));
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> sV(std::max(1, sieve->getMaxSupportSize()));
  ALE::ISieveVisitor::PointRetriever<SieveMesh::sieve_type> cV(std::max(1, sieve->getMaxConeSize()));
  for(PointSet::const_iterator fv_iter = fvBegin; fv_iter != fvEnd; ++fv_iter) {
    sieve->support(*fv_iter, sV);
    const SieveMesh::point_type *support = sV.getPoints();

    if (debug) std::cout << "Checking fault vertex " << *fv_iter << std::endl;
    const int sVsize = sV.getSize();
    for (int i=0; i < sVsize; ++i) {
      const int s = (!flipFault) ? i : sVsize - i - 1;
      sieve->cone(support[s], cV);
      const SieveMesh::point_type *cone = cV.getPoints();

      if (debug) std::cout << "  Checking cell " << support[s] << std::endl;
      if (faultCells.find(support[s]) != faultCells.end()) {
        cV.clear();
        continue;
      }
      face->clear();
      for(int c = 0; c < cV.getSize(); ++c) {
        if (faultVertices.find(cone[c]) != fvEnd) {
          if (debug) std::cout << "    contains fault vertex " << cone[c] << std::endl;
          face->insert(face->end(), cone[c]);
        } // if
      } // for
      if (face->size() > faceSize)
        throw ALE::Exception("Invalid fault mesh: Too many vertices of an "
                             "element on the fault");
      if (face->size() == faceSize) {
        if (debug) std::cout << "  Contains a face on the fault" << std::endl;
        ALE::Obj<SieveMesh::sieve_type::supportSet> preFace;
        if (dim < 2) {
          preFace = faultSieve->nJoin1(face);
        } else {
          preFace = faultSieve->nJoin(face, dim);
        }

        if (preFace->size() > 1) {
          throw ALE::Exception("Invalid fault sieve: Multiple faces from vertex set");
        } else if (preFace->size() == 1) {
          // Add the other cell neighbor for this face
          if (dim == 0) {
            faultSieve->addArrow(*faceVertices.begin(), support[s]);
          } else {
            faultSieve->addArrow(*preFace->begin(), support[s]);
          }
        } else if (preFace->size() == 0) {
          if (debug) std::cout << "  Orienting face " << f << std::endl;
          selection::getOrientedFace(mesh, support[s], face, numCorners, indices, &origVertices, &faceVertices);
          bdVertices[dim].clear();
          for(PointArray::const_iterator v_iter = faceVertices.begin(); v_iter != faceVertices.end(); ++v_iter) {
            bdVertices[dim].push_back(*v_iter);
            if (debug) std::cout << "    Boundary vertex " << *v_iter << std::endl;
          }
          if (dim == 0) {
            f = *faceVertices.begin();
          }

	  //std::cout << "dim: " << dim << ", faceSize: " << faceSize << ", numCorners: " << numCorners << std::endl;

          if (2 == dim && 4 == faceSize){
            if (debug) std::cout << "  Adding hex face " << f << std::endl;
            ALE::SieveBuilder<FlexMesh>::buildHexFaces(
		     faultSieve, orientation, dim, curElement, 
		     bdVertices, oFaultFaces, f, o);
          } else if ((1 == dim && 3 == faceSize) ||
		     (2 == dim && 9 == faceSize)){
            if (debug) std::cout << "  Adding quadratic hex face " << f
				 << std::endl;
            ALE::SieveBuilder<FlexMesh>::buildQuadraticHexFaces(
		     faultSieve, orientation, dim, curElement, 
		     bdVertices, oFaultFaces, f, o);
          } else if ((1 == dim && 3 == faceSize) ||
		     (2 == dim && 6 == faceSize)){
            if (debug) std::cout << "  Adding quadratic tri face " << f
				 << std::endl;
            ALE::SieveBuilder<FlexMesh>::buildQuadraticTetFaces(
		     faultSieve, orientation, dim, curElement, 
		     bdVertices, oFaultFaces, f, o);
          } else {
            if (debug) std::cout << "  Adding simplicial face " << f << std::endl;
            ALE::SieveBuilder<FlexMesh>::buildFaces(
		     faultSieve, orientation, dim, curElement,
		     bdVertices, oFaultFaces, f, o);
          }
          faultSieve->addArrow(f, support[s]);
          //faultSieve->view("");
          f++;
        } // if/else
        faultCells.insert(support[s]);
      } // if
      cV.clear();
    } // for
    sV.clear();
  } // for
  if (!faultSieve->commRank()) delete [] indices;
}

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::createFaultSieveFromFaces(const int dim,
                                                            const int firstCell,
                                                            const int numFaces,
                                                            const int faultVertices[],
                                                            const int faultCells[],
                                                            const ALE::Obj<SieveMesh>& mesh,
                                                            const ALE::Obj<FlexMesh::arrow_section_type>& orientation,
                                                            const ALE::Obj<FlexMesh::sieve_type>& faultSieve)
{
  typedef ALE::Selection<FlexMesh> selection;
  int                       faceSize   = 0; // The number of vertices in a mesh face
  int                       curCell    = firstCell;
  int                       curVertex  = 0;
  int                       newElement = curCell + dim*numFaces;
  int                       o          = 1;
  int                       f          = firstCell;
  const int                 debug      = mesh->debug();
  std::map<int,int*>        curElement;
  std::map<int,PointArray>  bdVertices;
  std::map<int,oPointArray> oFaultFaces;

  if (!faultSieve->commRank()) {
    faceSize = selection::numFaceVertices(mesh);
  }

  curElement[0]   = &curVertex;
  curElement[dim] = &curCell;
  for(int d = 1; d < dim; d++) {
    curElement[d] = &newElement;
  }

  // Loop over fault faces
  for(int face = 0; face < numFaces; ++face) {
    // Push oriented vertices of face
    bdVertices[dim].clear();
    for(int i = 0; i < faceSize; ++i) {
      bdVertices[dim].push_back(faultVertices[face*faceSize+i]);
      if (debug) std::cout << "    Boundary vertex " << faultVertices[face*faceSize+i] << std::endl;
    }
    // Create face
    if (faceSize != dim+1) {
      if (debug) std::cout << "  Adding hex face " << f << std::endl;
      ALE::SieveBuilder<FlexMesh>::buildHexFaces(faultSieve, orientation, dim, curElement, bdVertices, oFaultFaces, f, o);
    } else {
      if (debug) std::cout << "  Adding simplicial face " << f << std::endl;
      ALE::SieveBuilder<FlexMesh>::buildFaces(faultSieve, orientation, dim, curElement, bdVertices, oFaultFaces, f, o);
    }
    // Add arrow to cells
    faultSieve->addArrow(face, faultCells[face*2+0]);
    faultSieve->addArrow(face, faultCells[face*2+1]);
  }
}

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::orientFaultSieve(const int dim,
                                                   const ALE::Obj<SieveMesh>& mesh,
                                                   const ALE::Obj<FlexMesh::arrow_section_type>& orientation,
                                                   const ALE::Obj<FlexMesh>& fault)
{
  assert(!mesh.isNull());
  assert(!orientation.isNull());
  assert(!fault.isNull());

  typedef ALE::Selection<FlexMesh> selection;

  // Must check the orientation here
  const ALE::Obj<FlexMesh::sieve_type>& faultSieve = fault->getSieve();
  assert(!faultSieve.isNull());
  const SieveMesh::point_type firstFaultCell  = 
    *fault->heightStratum(1)->begin();
  const ALE::Obj<FlexMesh::label_sequence>& fFaces = fault->heightStratum(2);
  assert(!fFaces.isNull());
  const int numFaultFaces = fFaces->size();
  const int faultDepth = fault->depth()-1; // Depth of fault cells
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int faultFaceSize = 0; // The number of vertices in a face between fault cells
  int faceSize = 0; // The number of vertices in a mesh face
  const int debug = fault->debug();
  ALE::Obj<PointSet> newCells = new PointSet();
  assert(!newCells.isNull());
  ALE::Obj<PointSet> loopCells = new PointSet();
  assert(!loopCells.isNull());
  PointSet flippedCells;   // Incorrectly oriented fault cells
  PointSet facesSeen;      // Fault faces already considered
  PointSet cellsSeen;      // Fault cells already matched
  PointArray faceVertices;

  if (!fault->commRank()) {
    faceSize        = selection::numFaceVertices(mesh);
    numFaultCorners = faultSieve->nCone(firstFaultCell, faultDepth)->size();
    if (debug) std::cout << "  Fault corners " << numFaultCorners << std::endl;
    if (dim > 0) {
      assert(numFaultCorners == faceSize);
    } else {
      // dim is 0
      assert(numFaultCorners == faceSize-1);
    }
    if (faultDepth == 1) {
      faultFaceSize = 1;
    } else {
      if (dim > 0) {
	assert(fFaces->size() > 0);
	assert(faultDepth > 0);
	faultFaceSize = faultSieve->nCone(*fFaces->begin(), faultDepth-1)->size();
      } else {
	// dim is 0
	faultFaceSize = 1;
      } // if/else
    }
  }
  if (debug) std::cout << "  Fault face size " << faultFaceSize << std::endl;

  newCells->insert(firstFaultCell);
  while(facesSeen.size() != numFaultFaces) {
    ALE::Obj<PointSet> tmp = newCells; newCells = loopCells; loopCells = tmp;
        
    newCells->clear();
    if (!loopCells->size()) {throw ALE::Exception("Fault surface not a single connected component.");}
    // Loop over new cells
    for(PointSet::iterator c_iter = loopCells->begin(); c_iter != loopCells->end(); ++c_iter) {
      // Loop over edges of this cell
      const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>&     cone   = faultSieve->cone(*c_iter);
      const FlexMesh::sieve_type::traits::coneSequence::iterator eBegin = cone->begin();
      const FlexMesh::sieve_type::traits::coneSequence::iterator eEnd   = cone->end();

      for(FlexMesh::sieve_type::traits::coneSequence::iterator e_iter = eBegin; e_iter != eEnd; ++e_iter) {
        if (facesSeen.find(*e_iter) != facesSeen.end()) continue;
        facesSeen.insert(*e_iter);
        if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
        const ALE::Obj<FlexMesh::sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
        FlexMesh::sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

        // Throw out boundary fault faces
        if (support->size() < 2) continue;
        FlexMesh::point_type cellA = *s_iter; ++s_iter;
        FlexMesh::point_type cellB = *s_iter;
        bool flippedA = (flippedCells.find(cellA) != flippedCells.end());
        bool flippedB = (flippedCells.find(cellB) != flippedCells.end());
        bool seenA    = (cellsSeen.find(cellA) != cellsSeen.end());
        bool seenB    = (cellsSeen.find(cellB) != cellsSeen.end());

        if (!seenA) newCells->insert(cellA);
        if (!seenB) newCells->insert(cellB);
        if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
        // In 1D, just check that vertices match
        if (dim == 1) {
          const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
          FlexMesh::sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
          const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
          FlexMesh::sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
          int posA, posB;

          for(posA = 0; posA < 2; ++posA, ++iterA) if (*iterA == *e_iter) break;
          for(posB = 0; posB < 2; ++posB, ++iterB) if (*iterB == *e_iter) break;
          if (debug) std::cout << "    with face positions " << posA << " and " << posB << std::endl;
          if ((posA == 2) || (posB == 2)) {throw ALE::Exception("Could not find fault face in cone");}
          if ((posA == posB) ^ (flippedA || flippedB)) {
            if (debug) {
              std::cout << "Invalid orientation in fault mesh" << std::endl;
              std::cout << "  fault face: " << *e_iter << "  cellA: " << cellA << "  cellB: " << cellB << std::endl;
            }
            if (flippedA && flippedB) {throw ALE::Exception("Attempt to flip already flipped cell: Fault mesh is non-orientable");}
            if (seenA    && seenB)    {throw ALE::Exception("Previous cells do not match: Fault mesh is non-orientable");}
            if (!seenA && !flippedA) {
              flippedCells.insert(cellA);
            } else if (!seenB && !flippedB) {
              flippedCells.insert(cellB);
            } else {
              throw ALE::Exception("Inconsistent mesh orientation: Fault mesh is non-orientable");
            }
          }
        } else if (dim == 2) {
          // Check orientation
          ALE::MinimalArrow<FlexMesh::sieve_type::point_type,FlexMesh::sieve_type::point_type> arrowA(*e_iter, cellA);
          const int oA = orientation->restrictPoint(arrowA)[0];
          ALE::MinimalArrow<FlexMesh::sieve_type::point_type,FlexMesh::sieve_type::point_type> arrowB(*e_iter, cellB);
          const int oB = orientation->restrictPoint(arrowB)[0];
          const bool mismatch = (oA == oB);

          // Truth Table
          // mismatch    flips   action   mismatch   flipA ^ flipB   action
          //    F       0 flips    no        F             F           F
          //    F       1 flip     yes       F             T           T
          //    F       2 flips    no        T             F           T
          //    T       0 flips    yes       T             T           F
          //    T       1 flip     no
          //    T       2 flips    yes
          if (mismatch ^ (flippedA ^ flippedB)) {
            if (debug) {
              std::cout << "Invalid orientation in fault mesh" << std::endl;
              std::cout << "  fault face: " << *e_iter << "  cellA: " << cellA << "  cellB: " << cellB << std::endl;
            }
            if (flippedA && flippedB) {throw ALE::Exception("Attempt to flip already flipped cell: Fault mesh is non-orientable");}
            if (seenA    && seenB)    {throw ALE::Exception("Previous cells do not match: Fault mesh is non-orientable");}
            if (!seenA && !flippedA) {
              flippedCells.insert(cellA);
              if (debug) {std::cout << "    Scheduling cell " << cellA << " for flipping" << std::endl;}
            } else if (!seenB && !flippedB) {
              flippedCells.insert(cellB);
              if (debug) {std::cout << "    Scheduling cell " << cellB << " for flipping" << std::endl;}
            } else {
              throw ALE::Exception("Inconsistent mesh orientation: Fault mesh is non-orientable");
            }
          }
        }
        cellsSeen.insert(cellA);
        cellsSeen.insert(cellB);
      }
    }
  }
  for(PointSet::const_iterator f_iter = flippedCells.begin(); f_iter != flippedCells.end(); ++f_iter) {
    if (debug) std::cout << "  Reversing fault face " << *f_iter << std::endl;
    faceVertices.clear();
    const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>& cone = faultSieve->cone(*f_iter);
    for(FlexMesh::sieve_type::traits::coneSequence::iterator v_iter = cone->begin(); v_iter != cone->end(); ++v_iter) {
      faceVertices.insert(faceVertices.begin(), *v_iter);
    }
    faultSieve->clearCone(*f_iter);
    int color = 0;
    for(PointArray::const_iterator v_iter = faceVertices.begin(); v_iter != faceVertices.end(); ++v_iter) {
      faultSieve->addArrow(*v_iter, *f_iter, color++);
    }

    if (dim > 1) {
      // Here, they are edges, not vertices
      for(PointArray::const_iterator e_iter = faceVertices.begin(); e_iter != faceVertices.end(); ++e_iter) {
        ALE::MinimalArrow<FlexMesh::sieve_type::point_type,FlexMesh::sieve_type::point_type> arrow(*e_iter, *f_iter);
        int o = orientation->restrictPoint(arrow)[0];

        if (debug) std::cout << "    Reversing orientation of " << *e_iter <<"-->"<<*f_iter << " from " << o << " to " << -(o+1) << std::endl;
        o = -(o+1);
        orientation->updatePoint(arrow, &o);
      }
    }
  }
  flippedCells.clear();
  const FlexMesh::label_sequence::iterator fFacesBegin = fFaces->begin();
  const FlexMesh::label_sequence::iterator fFacesEnd = fFaces->end();
  for(FlexMesh::label_sequence::iterator e_iter = fFacesBegin; e_iter != fFacesEnd; ++e_iter) {
    if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
    // for each face get the support (2 fault cells)
    const ALE::Obj<FlexMesh::sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
    FlexMesh::sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

    // Throw out boundary fault faces
    if (support->size() > 1) {
      FlexMesh::point_type cellA = *s_iter; ++s_iter;
      FlexMesh::point_type cellB = *s_iter;

      if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
      // In 1D, just check that vertices match
      if (dim == 1) {
        const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
        FlexMesh::sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
        const ALE::Obj<FlexMesh::sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
        FlexMesh::sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
        int posA, posB;

        for(posA = 0; posA < 2; ++posA, ++iterA) if (*iterA == *e_iter) break;
        for(posB = 0; posB < 2; ++posB, ++iterB) if (*iterB == *e_iter) break;
        if (debug) std::cout << "    with face positions " << posA << " and " << posB << std::endl;
        if ((posA == 2) || (posB == 2)) {throw ALE::Exception("Could not find fault face in cone");}
        if (posA == posB) {
          std::cout << "Invalid orientation in fault mesh" << std::endl;
          std::cout << "  fault face: " << *e_iter << "  cellA: " << cellA << "  cellB: " << cellB << std::endl;
          throw ALE::Exception("Invalid orientation in fault mesh");
        }
      } else {
        // Check orientation
        ALE::MinimalArrow<FlexMesh::sieve_type::point_type,FlexMesh::sieve_type::point_type> arrowA(*e_iter, cellA);
        const int oA = orientation->restrictPoint(arrowA)[0];
        ALE::MinimalArrow<FlexMesh::sieve_type::point_type,FlexMesh::sieve_type::point_type> arrowB(*e_iter, cellB);
        const int oB = orientation->restrictPoint(arrowB)[0];

        if (oA == oB) {
          std::cout << "Invalid orientation in fault mesh" << std::endl;
          std::cout << "  fault face: " << *e_iter << "  cellA: " << cellA << "  cellB: " << cellB << std::endl;
          throw ALE::Exception("Invalid orientation in fault mesh");
        }
      }
    }
  }
  if (debug) fault->view("Oriented Fault mesh");
}


// End of file
