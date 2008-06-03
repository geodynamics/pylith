// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "CohesiveTopology.hh" // implementation of object methods
#include <Selection.hh> // Algorithms for submeshes

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::create(ALE::Obj<Mesh>* ifault,
					 const ALE::Obj<Mesh>& mesh,
					 const ALE::Obj<Mesh::int_section_type>& groupField,
					 const int materialId,
					 const bool constraintCell)
{ // create
  assert(0 != ifault);

  typedef ALE::SieveAlg<ALE::Mesh>  sieveAlg;
  typedef ALE::Selection<ALE::Mesh> selection;

  const int_section_type::chart_type& chart = groupField->getChart();
  PointSet faultVertices; // Vertices on fault
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  *ifault = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  const ALE::Obj<Mesh::sieve_type> ifaultSieve = new Mesh::sieve_type(sieve->comm(), 
                                                                      sieve->debug());
  ALE::Obj<ALE::Mesh> fault = new ALE::Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  ALE::Obj<ALE::Mesh::sieve_type> faultSieve = new ALE::Mesh::sieve_type(sieve->comm(), 
                                                                         sieve->debug());
  const int  depth      = mesh->depth();
  const int  numCells   = mesh->heightStratum(0)->size();
  int        numCorners = 0;    // The number of vertices in a mesh cell
  int        faceSize   = 0;    // The number of vertices in a mesh face
  int       *indices    = NULL; // The indices of a face vertex set in a cell
  int        oppositeVertex;    // For simplices, the vertex opposite a given face
  PointArray origVertices;
  PointArray faceVertices;
  PointArray neighborVertices;

  if (!fault->commRank()) {
    numCorners = mesh->getNumCellCorners();
    faceSize   = selection::numFaceVertices(mesh);
    indices    = new int[faceSize];
  }

  // Create set with vertices on fault
  for(int_section_type::chart_type::const_iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
    assert(!mesh->depth(*c_iter));
    if (groupField->getFiberDimension(*c_iter)) {faultVertices.insert(*c_iter);}
  } // for

  const PointSet::const_iterator fvBegin = faultVertices.begin();
  const PointSet::const_iterator fvEnd   = faultVertices.end();

  int f     = sieve->getBaseSize() + sieve->getCapSize();
  int debug = mesh->debug();
  ALE::Obj<PointSet> face = new PointSet();
  PointSet faultCells;
  
  // Create a sieve which captures the fault
  const int fDim = fault->getDimension();
  const Obj<Mesh::arrow_section_type>& orientation = fault->getArrowSection("orientation");
  std::map<int,int*>        curElement;
  std::map<int,PointArray>  bdVertices;
  std::map<int,PointArray>  faultFaces;
  std::map<int,oPointArray> oFaultFaces;
  int                       curCell    = f;
  int                       curVertex  = 0;
  int                       newElement = curCell + fDim*faultVertices.size();
  int                       o          = 1;

  curElement[0]   = &curVertex;
  curElement[fDim] = &curCell;
  for(int d = 1; d < fDim; d++) {
    curElement[d] = &newElement;
  }

  // This only works for uninterpolated meshes
  assert((mesh->depth() == 1) || (mesh->depth() == -1));
  ALE::ISieveVisitor::PointRetriever<sieve_type> sV(std::max(1, sieve->getMaxSupportSize()));
  ALE::ISieveVisitor::PointRetriever<sieve_type> cV(std::max(1, sieve->getMaxConeSize()));
  for(PointSet::const_iterator fv_iter = fvBegin; fv_iter != fvEnd; ++fv_iter) {
    sieve->support(*fv_iter, sV);
    const Mesh::point_type *support = sV.getPoints();

    if (debug) std::cout << "Checking fault vertex " << *fv_iter << std::endl;
    for(int s = 0; s < sV.getSize(); ++s) {
      sieve->cone(support[s], cV);
      const Mesh::point_type *cone = cV.getPoints();

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
        ALE::Obj<sieve_type::supportSet> preFace;
        if (fDim < 2) {
          preFace = faultSieve->nJoin1(face);
        } else {
          preFace = faultSieve->nJoin(face, fDim);
        }

        if (preFace->size() > 1) {
          throw ALE::Exception("Invalid fault sieve: Multiple faces from vertex set");
        } else if (preFace->size() == 1) {
          // Add the other cell neighbor for this face
          if (fDim == 0) {
            faultSieve->addArrow(*faceVertices.begin(), support[s]);
          } else {
            faultSieve->addArrow(*preFace->begin(), support[s]);
          }
        } else if (preFace->size() == 0) {
          if (debug) std::cout << "  Orienting face " << f << std::endl;
          selection::getOrientedFace(mesh, support[s], face, numCorners, indices, &origVertices, &faceVertices);
          bdVertices[fDim].clear();
          for(PointArray::const_iterator v_iter = faceVertices.begin(); v_iter != faceVertices.end(); ++v_iter) {
            bdVertices[fDim].push_back(*v_iter);
            if (debug) std::cout << "    Boundary vertex " << *v_iter << std::endl;
          }
          if (fDim == 0) {
            f = *faceVertices.begin();
          }
          if (faceSize != fDim+1) {
            if (debug) std::cout << "  Adding hex face " << f << std::endl;
            ALE::SieveBuilder<ALE::Mesh>::buildHexFaces(faultSieve, orientation, fDim, curElement, bdVertices, oFaultFaces, f, o);
          } else {
            if (debug) std::cout << "  Adding simplicial face " << f << std::endl;
            ALE::SieveBuilder<ALE::Mesh>::buildFaces(faultSieve, orientation, fDim, curElement, bdVertices, oFaultFaces, f, o);
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
  fault->setSieve(faultSieve);
  fault->stratify();
  faultCells.clear();
  if (debug) fault->view("Fault mesh");
  const ALE::Obj<ALE::Mesh> faultBd = ALE::Selection<ALE::Mesh>::boundary(fault);
  if (debug) faultBd->view("Fault boundary mesh");

  // Orient the fault sieve
  // Must check the orientation here
  const Mesh::point_type firstFaultCell = *fault->heightStratum(1)->begin();
  const ALE::Obj<Mesh::label_sequence>& fFaces = fault->heightStratum(2);
  const int              numFaultFaces  = fFaces->size();
  int faultDepth      = fault->depth()-1; // Depth of fault cells
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int faultFaceSize   = 0; // The number of vertices in a face between fault cells
  PointSet flippedCells;   // Incorrectly oriented fault cells
  PointSet facesSeen;      // Fault faces already considered
  PointSet cellsSeen;      // Fault cells already matched
  Obj<PointSet> newCells  = new PointSet();
  Obj<PointSet> loopCells = new PointSet();

  if (!fault->commRank()) {
    numFaultCorners = faultSieve->nCone(firstFaultCell, faultDepth)->size();
    if (debug) std::cout << "  Fault corners " << numFaultCorners << std::endl;
    if (fDim == 0) {
      assert(numFaultCorners == faceSize-1);
    } else {
      assert(numFaultCorners == faceSize);
    }
    if (faultDepth == 1) {
      faultFaceSize = 1;
    } else {
      faultFaceSize = faultSieve->nCone(*fFaces->begin(), faultDepth-1)->size();
    }
  }
  if (debug) std::cout << "  Fault face size " << faultFaceSize << std::endl;

  newCells->insert(firstFaultCell);
  while(facesSeen.size() != numFaultFaces) {
    Obj<PointSet> tmp = newCells; newCells = loopCells; loopCells = tmp;
        
    newCells->clear();
    if (!loopCells->size()) {throw ALE::Exception("Fault surface not a single connected component.");}
    // Loop over new cells
    for(PointSet::iterator c_iter = loopCells->begin(); c_iter != loopCells->end(); ++c_iter) {
      // Loop over edges of this cell
      const Obj<ALE::Mesh::sieve_type::traits::coneSequence>&     cone   = faultSieve->cone(*c_iter);
      const ALE::Mesh::sieve_type::traits::coneSequence::iterator eBegin = cone->begin();
      const ALE::Mesh::sieve_type::traits::coneSequence::iterator eEnd   = cone->end();

      for(ALE::Mesh::sieve_type::traits::coneSequence::iterator e_iter = eBegin; e_iter != eEnd; ++e_iter) {
        if (facesSeen.find(*e_iter) != facesSeen.end()) continue;
        facesSeen.insert(*e_iter);
        if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
        const Obj<ALE::Mesh::sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
        ALE::Mesh::sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

        // Throw out boundary fault faces
        if (support->size() < 2) continue;
        ALE::Mesh::point_type cellA = *s_iter; ++s_iter;
        ALE::Mesh::point_type cellB = *s_iter;
        bool flippedA = (flippedCells.find(cellA) != flippedCells.end());
        bool flippedB = (flippedCells.find(cellB) != flippedCells.end());
        bool seenA    = (cellsSeen.find(cellA) != cellsSeen.end());
        bool seenB    = (cellsSeen.find(cellB) != cellsSeen.end());

        if (!seenA) newCells->insert(cellA);
        if (!seenB) newCells->insert(cellB);
        if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
        // In 1D, just check that vertices match
        if (fDim == 1) {
          const Obj<ALE::Mesh::sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
          ALE::Mesh::sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
          const Obj<ALE::Mesh::sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
          ALE::Mesh::sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
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
        } else if (fDim == 2) {
          // Check orientation
          ALE::MinimalArrow<ALE::Mesh::sieve_type::point_type,ALE::Mesh::sieve_type::point_type> arrowA(*e_iter, cellA);
          const int oA = orientation->restrictPoint(arrowA)[0];
          ALE::MinimalArrow<ALE::Mesh::sieve_type::point_type,ALE::Mesh::sieve_type::point_type> arrowB(*e_iter, cellB);
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
    const ALE::Obj<ALE::Mesh::sieve_type::traits::coneSequence>& cone = faultSieve->cone(*f_iter);
    for(ALE::Mesh::sieve_type::traits::coneSequence::iterator v_iter = cone->begin(); v_iter != cone->end(); ++v_iter) {
      faceVertices.insert(faceVertices.begin(), *v_iter);
    }
    faultSieve->clearCone(*f_iter);
    int color = 0;
    for(PointArray::const_iterator v_iter = faceVertices.begin(); v_iter != faceVertices.end(); ++v_iter) {
      faultSieve->addArrow(*v_iter, *f_iter, color++);
    }

    if (fDim > 1) {
      // Here, they are edges, not vertices
      for(PointArray::const_iterator e_iter = faceVertices.begin(); e_iter != faceVertices.end(); ++e_iter) {
        ALE::MinimalArrow<ALE::Mesh::sieve_type::point_type,ALE::Mesh::sieve_type::point_type> arrow(*e_iter, *f_iter);
        int o = orientation->restrictPoint(arrow)[0];

        if (debug) std::cout << "    Reversing orientation of " << *e_iter <<"-->"<<*f_iter << " from " << o << " to " << -(o+1) << std::endl;
        o = -(o+1);
        orientation->updatePoint(arrow, &o);
      }
    }
  }
  flippedCells.clear();
  for(ALE::Mesh::label_sequence::iterator e_iter = fFaces->begin(); e_iter != fFaces->end(); ++e_iter) {
    if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
    // for each face get the support (2 fault cells)
    const Obj<ALE::Mesh::sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
    ALE::Mesh::sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

    // Throw out boundary fault faces
    if (support->size() > 1) {
      ALE::Mesh::point_type cellA = *s_iter; ++s_iter;
      ALE::Mesh::point_type cellB = *s_iter;

      if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
      // In 1D, just check that vertices match
      if (fDim == 1) {
        const Obj<ALE::Mesh::sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
        ALE::Mesh::sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
        const Obj<ALE::Mesh::sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
        ALE::Mesh::sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
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
        ALE::MinimalArrow<ALE::Mesh::sieve_type::point_type,ALE::Mesh::sieve_type::point_type> arrowA(*e_iter, cellA);
        const int oA = orientation->restrictPoint(arrowA)[0];
        ALE::MinimalArrow<ALE::Mesh::sieve_type::point_type,ALE::Mesh::sieve_type::point_type> arrowB(*e_iter, cellB);
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

  // Convert fault to an IMesh
  Mesh::renumbering_type& renumbering = (*ifault)->getRenumbering();
  (*ifault)->setSieve(ifaultSieve);
  ALE::ISieveConverter::convertMesh(*fault, *(*ifault), renumbering, false);
  renumbering.clear();
  fault      = NULL;
  faultSieve = NULL;

  // Add new shadow vertices and possibly Lagrange multipler vertices
  const ALE::Obj<Mesh::label_sequence>& fVertices = (*ifault)->depthStratum(0);
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<std::set<std::string> >& groupNames = mesh->getIntSections();
  Mesh::point_type newPoint = sieve->getBaseSize() + sieve->getCapSize();
  const int        numFaultVertices = fVertices->size();
  std::map<Mesh::point_type,Mesh::point_type> vertexRenumber;
  std::map<Mesh::point_type,Mesh::point_type> cellRenumber;

  // I know how to fix
  for(Mesh::label_sequence::iterator v_iter = fVertices->begin(); v_iter != fVertices->end(); ++v_iter, ++newPoint) {
    vertexRenumber[*v_iter] = newPoint;
    if (debug) std::cout << "Duplicating " << *v_iter << " to " << vertexRenumber[*v_iter] << std::endl;

    // Add shadow and constraint vertices (if they exist) to group
    // associated with fault
    groupField->addPoint(newPoint, 1);
    if (constraintCell)
      groupField->addPoint(newPoint+numFaultVertices, 1);

    // Add shadow vertices to other groups, don't add constraint
    // vertices (if they exist) because we don't want BC, etc to act
    // on constraint vertices
    for(std::set<std::string>::const_iterator name = groupNames->begin();
       name != groupNames->end(); ++name) {
      const ALE::Obj<int_section_type>& group = mesh->getIntSection(*name);
      if (group->getFiberDimension(*v_iter))
        group->addPoint(newPoint, 1);
    } // for
  } // for
  for(std::set<std::string>::const_iterator name = groupNames->begin();
      name != groupNames->end(); ++name) {
    mesh->reallocate(mesh->getIntSection(*name));
  } // for
  if (constraintCell) newPoint += numFaultVertices;

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<Mesh::label_sequence>& faces = (*ifault)->heightStratum(1);
  const ALE::Obj<Mesh::label_type>& material = mesh->getLabel("material-id");
  const int firstCohesiveCell = newPoint;
  PointSet replaceCells;
  PointSet noReplaceCells;
  PointSet replaceVertices;
  ALE::ISieveVisitor::PointRetriever<sieve_type> sV2(std::max(1, ifaultSieve->getMaxSupportSize()));
  ALE::ISieveVisitor::NConeRetriever<sieve_type> cV2(*ifaultSieve, (size_t) pow(std::max(1, ifaultSieve->getMaxConeSize()), (*ifault)->depth()));

  for(Mesh::label_sequence::iterator f_iter = faces->begin(); f_iter != faces->end(); ++f_iter, ++newPoint) {
    const Mesh::point_type face = *f_iter;
    if (debug) std::cout << "Considering fault face " << face << std::endl;
    ifaultSieve->support(face, sV2);
    const Mesh::point_type *cells = sV2.getPoints();
    Mesh::point_type cell = cells[0];
    Mesh::point_type otherCell;

    if (debug) std::cout << "  Checking orientation against cell " << cell << std::endl;
    selection::getOrientedFace(mesh, cell, &vertexRenumber, numCorners, indices, &origVertices, &faceVertices);

    ALE::ISieveTraversal<sieve_type>::orientedClosure(*ifaultSieve, face, cV2);
    const int               coneSize = cV2.getSize();
    const Mesh::point_type *faceCone = cV2.getPoints();
    //ifaultSieve->cone(face, cV2);
    //const int               coneSize = cV2.getSize() ? cV2.getSize()   : 1;
    //const Mesh::point_type *faceCone = cV2.getSize() ? cV2.getPoints() : &face;
    bool                    found    = true;

    if (numFaultCorners == 0) {
      found = false;
    } else if (numFaultCorners == 2) {
      if (faceVertices[0] != faceCone[0]) found = false;
    } else {
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVertices[v] != faceCone[0])) ++v;
      for(int c = 0; c < coneSize; ++c, ++v) {
        if (debug) std::cout << "    Checking " << faceCone[c] << " against " << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != faceCone[c]) {
          found = false;
          break;
        }
      }
    }
    if (found) {
      if (debug) std::cout << "  Choosing other cell" << std::endl;
      otherCell = cell;
      cell = cells[1];
    } else {
      otherCell = cells[1];
      if (debug) std::cout << "  Verifing reverse orientation" << std::endl;
      found = true;
      int v = 0;
      if (numFaultCorners > 0) {
        // Locate first vertex
        while((v < numFaultCorners) && (faceVertices[v] != faceCone[coneSize-1])) ++v;
        for(int c = coneSize-1; c >= 0; --c, ++v) {
          if (debug) std::cout << "    Checking " << faceCone[c] << " against " << faceVertices[v%numFaultCorners] << std::endl;
          if (faceVertices[v%numFaultCorners] != faceCone[c]) {
            found = false;
            break;
          }
        }
      }
      assert(found);
    }
    noReplaceCells.insert(otherCell);
    replaceCells.insert(cell);
    replaceVertices.insert(faceCone, &faceCone[coneSize]);
    cellRenumber[cell] = newPoint;
    // Adding cohesive cell (not interpolated)
	if (debug) std::cout << "  Creating cohesive cell " << newPoint << std::endl;
    for(int c = 0; c < coneSize; ++c) {
      if (debug) std::cout << "    vertex " << faceCone[c] << std::endl;
      sieve->addArrow(faceCone[c], newPoint);
    }
    for(int c = 0; c < coneSize; ++c) {
      if (debug) std::cout << "    shadow vertex " << vertexRenumber[faceCone[c]] << std::endl;
      sieve->addArrow(vertexRenumber[faceCone[c]], newPoint);
    }
    if (constraintCell) {
      for(int c = 0; c < coneSize; ++c) {
        if (debug) std::cout << "    Lagrange vertex " << vertexRenumber[faceCone[c]]+numFaultVertices << std::endl;
        sieve->addArrow(vertexRenumber[faceCone[c]]+numFaultVertices, newPoint);
      }
    }
    mesh->setValue(material, newPoint, materialId);
    sV2.clear();
    cV2.clear();
  } // for
  // Add new arrows for support of replaced vertices
  for(PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != replaceVertices.end(); ++v_iter) {
    sieve->support(*v_iter, sV);
    const Mesh::point_type *support = sV.getPoints();

    for(int s = 0; s < sV.getSize(); ++s) {
      if (replaceCells.find(support[s]) != replaceCells.end()) {
        sieve->addArrow(vertexRenumber[*v_iter], support[s]);
      }
    }
    sV.clear();
  }
  sieve->reallocate();
  // More checking
  PointSet replaceCellsBase(replaceCells);

  const ALE::Obj<Mesh::label_sequence>& faultBdVerts = faultBd->depthStratum(0);
  PointSet faultBdVertices;

  faultBdVertices.insert(faultBdVerts->begin(), faultBdVerts->end());
  for(PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != replaceVertices.end(); ++v_iter) {
    if (faultBdVertices.find(*v_iter) != faultBdVertices.end()) continue;
    classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, faultBdVertices, replaceCells, noReplaceCells, debug);
  }
  for(PointSet::const_iterator v_iter = faultBdVertices.begin(); v_iter != faultBdVertices.end(); ++v_iter) {
    classifyCells(sieve, *v_iter, depth, faceSize, firstCohesiveCell, faultBdVertices, replaceCells, noReplaceCells, debug);
  }

  // More checking
  const bool                         firstFault    = !mesh->hasRealSection("replacedCells");
  const ALE::Obj<real_section_type>& replacedCells = mesh->getRealSection("replacedCells");
  PointSet cellNeighbors;

  if (firstFault) {
    const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

    replacedCells->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()), *std::max_element(cells->begin(), cells->end())+1));
    replacedCells->setFiberDimension(cells, 1);
    replacedCells->allocatePoint();
  }
  for(PointSet::const_iterator c_iter = noReplaceCells.begin(); c_iter != noReplaceCells.end(); ++c_iter) {
    const double minusOne = -1.0;

    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &minusOne);
    } else {
      const double minusTwo = -2.0;

      replacedCells->updatePoint(*c_iter, &minusTwo);
    }
  }
  for(PointSet::const_iterator c_iter = replaceCells.begin(); c_iter != replaceCells.end(); ++c_iter) {
    if (replaceCellsBase.find(*c_iter) != replaceCellsBase.end()) {
      const double one = 1.0;

      if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
        replacedCells->updatePoint(*c_iter, &one);
      } else {
        const double two = 2.0;

        replacedCells->updatePoint(*c_iter, &two);
      }
      continue;
    }
    const double ten = 10.0;

    if (replacedCells->restrictPoint(*c_iter)[0] == 0.0) {
      replacedCells->updatePoint(*c_iter, &ten);
    } else {
        const double twenty = 20.0;

        replacedCells->updatePoint(*c_iter, &twenty);
    }
    // There should be a way to check for boundary elements
    if (mesh->getDimension() == 1) {
      if (cellNeighbors.size() > 2) {
        std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
        throw ALE::Exception("Invalid number of neighbors");
      }
    } else if (mesh->getDimension() == 2) {
      if (numCorners == 3) {
        if (cellNeighbors.size() > 3) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
	}
      } else if (numCorners == 4) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        }
      }
    } else if (mesh->getDimension() == 3) {
      if (numCorners == 4) {
        if (cellNeighbors.size() > 4) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        }
      } else if (numCorners == 8) {
        if (cellNeighbors.size() > 6) {
          std::cout << "Cell " << *c_iter << " has an invalid number of neighbors " << cellNeighbors.size() << std::endl;
          throw ALE::Exception("Invalid number of neighbors");
        }
      }
    }
  }
  ReplaceVisitor<sieve_type,std::map<Mesh::point_type,Mesh::point_type> > rVc(vertexRenumber, std::max(1, sieve->getMaxConeSize()), debug);

  for(PointSet::const_iterator c_iter = replaceCells.begin(); c_iter != replaceCells.end(); ++c_iter) {
    sieve->cone(*c_iter, rVc);
    if (rVc.mappedPoint()) {
      if (debug) std::cout << "  Replacing cell " << *c_iter << std::endl;
      sieve->setCone(rVc.getPoints(), *c_iter);
    }
    rVc.clear();
  }
  ReplaceVisitor<sieve_type,std::map<Mesh::point_type,Mesh::point_type> > rVs(cellRenumber, std::max(1, sieve->getMaxSupportSize()), debug);

  for(PointSet::const_iterator v_iter = replaceVertices.begin(); v_iter != replaceVertices.end(); ++v_iter) {
    sieve->support(*v_iter, rVs);
    if (rVs.mappedPoint()) {
      if (debug) std::cout << "  Replacing support " << *v_iter << std::endl;
      sieve->setSupport(*v_iter, rVs.getPoints());
    }
    rVs.clear();
  }
  if (!(*ifault)->commRank()) delete [] indices;
  mesh->stratify();
  const std::string labelName("censored depth");

  if (!mesh->hasLabel(labelName)) {
    const ALE::Obj<Mesh::label_type>& label = mesh->createLabel(labelName);

    _computeCensoredDepth(label, mesh->getSieve(), firstCohesiveCell-(constraintCell?numFaultVertices:0));
  } else {
    // Insert new shadow vertices into existing label
    const ALE::Obj<Mesh::label_type>& label = mesh->getLabel(labelName);

    for(std::map<int,int>::const_iterator v_iter = vertexRenumber.begin(); v_iter != vertexRenumber.end(); ++v_iter) {
      mesh->setValue(label, v_iter->second, 0);
    }
  }
  if (debug) mesh->view("Mesh with Cohesive Elements");
  mesh->view("Mesh with Cohesive Elements");
  mesh->getLabel("depth")->view("Depth");
  mesh->getLabel("censored depth")->view("Censored Depth");

  // Fix coordinates
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& fVertices2 = (*ifault)->depthStratum(0);

  coordinates->view("Coordinates without shadow vertices");
  for(Mesh::label_sequence::iterator v_iter = fVertices2->begin();
      v_iter != fVertices2->end();
      ++v_iter) {
    coordinates->addPoint(vertexRenumber[*v_iter],
			  coordinates->getFiberDimension(*v_iter));
    if (constraintCell) {
      coordinates->addPoint(vertexRenumber[*v_iter]+numFaultVertices,
			  coordinates->getFiberDimension(*v_iter));
    }
  } // for
  mesh->reallocate(coordinates);
  for(Mesh::label_sequence::iterator v_iter = fVertices2->begin();
      v_iter != fVertices2->end();
      ++v_iter) {
    coordinates->updatePoint(vertexRenumber[*v_iter], 
			     coordinates->restrictPoint(*v_iter));
    if (constraintCell) {
      coordinates->updatePoint(vertexRenumber[*v_iter]+numFaultVertices,
			     coordinates->restrictPoint(*v_iter));
    }
  }
  if (debug) coordinates->view("Coordinates with shadow vertices");
  coordinates->view("Coordinates with shadow vertices");
} // createCohesiveCells

// ----------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::CohesiveTopology::createParallel(
		ALE::Obj<Mesh>* ifault,
		std::map<Mesh::point_type, Mesh::point_type>* cohesiveToFault,
		const ALE::Obj<Mesh>& mesh,
		const int materialId,
		const bool constraintCell)
{
  assert(0 != ifault);
  assert(0 != cohesiveToFault);

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  *ifault = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  const ALE::Obj<sieve_type> ifaultSieve = new sieve_type(sieve->comm(), sieve->debug());
  ALE::Obj<ALE::Mesh> fault = new ALE::Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  ALE::Obj<ALE::Mesh::sieve_type> faultSieve = new ALE::Mesh::sieve_type(sieve->comm(), sieve->debug());
  cohesiveToFault->clear();

  const ALE::Obj<Mesh::label_sequence>& cohesiveCells = mesh->getLabelStratum("material-id", materialId);
  const Mesh::label_sequence::iterator cBegin = cohesiveCells->begin();
  const Mesh::label_sequence::iterator cEnd = cohesiveCells->end();
  const int sieveEnd = sieve->getBaseSize() + sieve->getCapSize();
  const int numFaces = cohesiveCells->size();
  int globalSieveEnd = 0;
  int globalFaceOffset = 0;

  MPI_Allreduce((void *) &sieveEnd, (void *) &globalSieveEnd, 1, MPI_INT, MPI_SUM, sieve->comm());
  MPI_Scan((void *) &numFaces, (void *) &globalFaceOffset, 1, MPI_INT, MPI_SUM, sieve->comm());
  int face = globalSieveEnd + globalFaceOffset - numFaces;

  ALE::ISieveVisitor::PointRetriever<sieve_type> cV(sieve->getMaxConeSize());

  for(Mesh::label_sequence::iterator c_iter = cBegin; c_iter != cEnd; ++c_iter) {
    sieve->cone(*c_iter, cV);
    const int               coneSize = cV.getSize();
    const Mesh::point_type *cone     = cV.getPoints();
    int                     color    = 0;

    if (!constraintCell) {
      const int faceSize = coneSize / 2;
      assert(0 == coneSize % faceSize);

      // Use first vertices (negative side of the fault) for fault mesh
      for(int i = 0; i < faceSize; ++i) {
        faultSieve->addArrow(cone[i], face, color++);
      }
    } else {
      const int faceSize = coneSize / 3;
      assert(0 == coneSize % faceSize);

      // Use last vertices (contraints) for fault mesh
      for(int i = 2*faceSize; i < 3*faceSize; ++i) {
        faultSieve->addArrow(cone[i], face, color++);
      }
    } // if/else
    (*cohesiveToFault)[*c_iter] = face;
    ++face;
    cV.clear();
  } // for
  fault->setSieve(faultSieve);
  fault->stratify();

  // Convert fault to an IMesh
  std::map<Mesh::point_type,Mesh::point_type> renumbering;
  (*ifault)->setSieve(ifaultSieve);
  ALE::ISieveConverter::convertMesh(*fault, *(*ifault), renumbering, false);
  renumbering.clear();
  fault      = NULL;
  faultSieve = NULL;

  const ALE::Obj<Mesh::label_sequence>& faultCells = (*ifault)->heightStratum(0);
  assert(!faultCells.isNull());
  for(Mesh::label_sequence::iterator c_iter = cBegin, f_iter=faultCells->begin(); c_iter != cEnd; ++c_iter, ++f_iter) {
    (*cohesiveToFault)[*c_iter] = *f_iter;
  }
    
#if 1
  (*ifault)->setRealSection("coordinates", mesh->getRealSection("coordinates"));
#else
  const ALE::Obj<Mesh::real_section_type>& coordinates  = mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::real_section_type>& fCoordinates = (*fault)->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>&    vertices     = (*fault)->depthStratum(0);
  const Mesh::label_sequence::iterator     vBegin       = vertices->begin();
  const Mesh::label_sequence::iterator     vEnd         = vertices->end();

  for(Mesh::label_sequence::iterator v_iter = vBegin;
      v_iter != vEnd;
      ++v_iter) {
    fCoordinates->setFiberDimension(*v_iter, coordinates->getFiberDimension(*v_iter));
  }
  (*fault)->allocate(fCoordinates);
  for(Mesh::label_sequence::iterator v_iter = vBegin;
      v_iter != vEnd;
      ++v_iter) {
    fCoordinates->updatePoint(*v_iter, coordinates->restrictPoint(*v_iter));
  }
#endif

  // Create the parallel overlap
  //   Can I figure this out in a nicer way?
  Obj<Mesh::send_overlap_type> sendParallelMeshOverlap = (*ifault)->getSendOverlap();
  Obj<Mesh::recv_overlap_type> recvParallelMeshOverlap = (*ifault)->getRecvOverlap();
  ALE::SetFromMap<std::map<Mesh::point_type,Mesh::point_type> > globalPoints((*ifault)->getRenumbering());

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, (*ifault)->getRenumbering(), sendParallelMeshOverlap, recvParallelMeshOverlap);
  (*ifault)->setCalculatedOverlap(true);
}

// ----------------------------------------------------------------------
void
pylith::faults::CohesiveTopology::classifyCells(const ALE::Obj<Mesh::sieve_type>& sieve,
                                                const Mesh::point_type& vertex,
                                                const int depth,
                                                const int faceSize,
                                                const Mesh::point_type& firstCohesiveCell,
                                                const PointSet& faultBdVertices,
                                                PointSet& replaceCells,
                                                PointSet& noReplaceCells,
                                                const int debug)
{
  // Replace all cells on a given side of the fault with a vertex on the fault
  ClassifyVisitor<Mesh::sieve_type> cV(*sieve, replaceCells, noReplaceCells, firstCohesiveCell, faceSize, debug);
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
template<class InputPoints>
bool
pylith::faults::CohesiveTopology::_compatibleOrientation(const ALE::Obj<Mesh>& mesh,
                                                         const Mesh::point_type& p,
                                                         const Mesh::point_type& q,
                                                         const int numFaultCorners,
                                                         const int faultFaceSize,
                                                         const int faultDepth,
                                                         const Obj<InputPoints>& points,
                                                         int indices[],
                                                         PointArray *origVertices,
                                                         PointArray *faceVertices,
                                                         PointArray *neighborVertices)
{
  typedef ALE::Selection<Mesh> selection;
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
pylith::faults::CohesiveTopology::_computeCensoredDepth(const ALE::Obj<Mesh::label_type>& depth,
                                                        const ALE::Obj<Mesh::sieve_type>& sieve,
                                                        const Mesh::point_type& firstCohesiveCell)
{
  Mesh::DepthVisitor d(*sieve, firstCohesiveCell, *depth);

  sieve->roots(d);
  while(d.isModified()) {
    // FIX: Avoid the copy here somehow by fixing the traversal
    std::vector<Mesh::point_type> modifiedPoints(d.getModifiedPoints().begin(), d.getModifiedPoints().end());

    d.clear();
    sieve->support(modifiedPoints, d);
  }
};
// End of file
