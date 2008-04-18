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
pylith::faults::CohesiveTopology::create(ALE::Obj<Mesh>* fault,
					 const ALE::Obj<Mesh>& mesh,
					 const ALE::Obj<Mesh::int_section_type>& groupField,
					 const int materialId,
					 const bool constraintCell)
{ // create
  assert(0 != fault);

  typedef ALE::SieveAlg<Mesh>                    sieveAlg;
  typedef ALE::Selection<Mesh>                   selection;

  const int_section_type::chart_type& chart = groupField->getChart();
  PointSet faultVertices; // Vertices on fault
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  *fault = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  const ALE::Obj<sieve_type> faultSieve = new sieve_type(sieve->comm(), 
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

  if (!(*fault)->commRank()) {
    numCorners = sieve->nCone(*mesh->heightStratum(0)->begin(), depth)->size();
#if 1
    faceSize   = selection::numFaceVertices(mesh);
#else
    faceSize   = selection::numFaceVertices(*mesh->heightStratum(0)->begin(), mesh);
#endif
    indices    = new int[faceSize];
  }

  // Create set with vertices on fault
  for(int_section_type::chart_type::iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    assert(!mesh->depth(*c_iter));
    faultVertices.insert(*c_iter);
  } // for

  const PointSet::const_iterator fvBegin = faultVertices.begin();
  const PointSet::const_iterator fvEnd   = faultVertices.end();

  int f = sieve->base()->size() + sieve->cap()->size();
  int debug = mesh->debug();
  ALE::Obj<PointSet> face = new PointSet();
  PointSet faultCells;
  
  // Create a sieve which captures the fault
  const int fDim = (*fault)->getDimension();
  const Obj<Mesh::arrow_section_type>& orientation = (*fault)->getArrowSection("orientation");
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

  for(PointSet::const_iterator fv_iter = fvBegin; fv_iter != fvEnd; ++fv_iter) {
    const ALE::Obj<sieveAlg::supportArray>& cells =
      sieveAlg::nSupport(mesh, *fv_iter, depth);
    const sieveAlg::supportArray::iterator cBegin = cells->begin();
    const sieveAlg::supportArray::iterator cEnd   = cells->end();
    
    if (debug)
      std::cout << "Checking fault vertex " << *fv_iter << std::endl;
    for(sieveAlg::supportArray::iterator c_iter = cBegin;
        c_iter != cEnd;
        ++c_iter) {
      if (debug) std::cout << "  Checking cell " << *c_iter << std::endl;
      if (faultCells.find(*c_iter) != faultCells.end())	continue;
      const ALE::Obj<sieveAlg::coneArray>& cone =
        sieveAlg::nCone(mesh, *c_iter, mesh->height());
      const sieveAlg::coneArray::iterator vBegin = cone->begin();
      const sieveAlg::coneArray::iterator vEnd   = cone->end();

      face->clear();
      for(sieveAlg::coneArray::iterator v_iter = vBegin;
          v_iter != vEnd;
          ++v_iter) {
        if (faultVertices.find(*v_iter) != fvEnd) {
          if (debug) std::cout << "    contains fault vertex " << *v_iter << std::endl;
          face->insert(face->end(), *v_iter);
        } // if
      } // for
      if (face->size() > faceSize)
        throw ALE::Exception("Invalid fault mesh: Too many vertices of an "
                             "element on the fault");
      if (face->size() == faceSize) {
        if (debug)
          std::cout << "  Contains a face on the fault" << std::endl;
        ALE::Obj<sieve_type::supportSet> preFace;
        if (fDim < 2) {
          preFace = faultSieve->nJoin1(face);
        } else {
          preFace = faultSieve->nJoin(face, fDim);
        }

        if (preFace->size() > 1) {
          throw ALE::Exception("Invalid fault sieve: Multiple faces from "
                               "vertex set");
        } else if (preFace->size() == 1) {
          // Add the other cell neighbor for this face
          if (fDim == 0) {
            faultSieve->addArrow(*faceVertices.begin(), *c_iter);
          } else {
            faultSieve->addArrow(*preFace->begin(), *c_iter);
          }
        } else if (preFace->size() == 0) {
          if (debug) std::cout << "  Orienting face " << f << std::endl;
          selection::getOrientedFace(mesh, *c_iter, face, numCorners, indices, &origVertices, &faceVertices);
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
            ALE::SieveBuilder<Mesh>::buildHexFaces(faultSieve, orientation, fDim, curElement, bdVertices, oFaultFaces, f, o);
          } else {
            if (debug) std::cout << "  Adding simplicial face " << f << std::endl;
            ALE::SieveBuilder<Mesh>::buildFaces(faultSieve, orientation, fDim, curElement, bdVertices, oFaultFaces, f, o);
          }
          faultSieve->addArrow(f, *c_iter);
          //faultSieve->view("");
          f++;
        } // if/else
        faultCells.insert(*c_iter);
      } // if
    } // for
  } // for
  (*fault)->setSieve(faultSieve);
  (*fault)->stratify();
  faultCells.clear();
  if (debug) (*fault)->view("Fault mesh");
  const ALE::Obj<Mesh> faultBd = ALE::Selection<Mesh>::boundary(*fault);
  if (debug) faultBd->view("Fault boundary mesh");

  // Orient the fault sieve
  // Must check the orientation here
  const Mesh::point_type firstFaultCell = *(*fault)->heightStratum(1)->begin();
  const ALE::Obj<Mesh::label_sequence>& fFaces = (*fault)->heightStratum(2);
  const int              numFaultFaces  = fFaces->size();
  int faultDepth      = (*fault)->depth()-1; // Depth of fault cells
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int faultFaceSize   = 0; // The number of vertices in a face between fault cells
  PointSet flippedCells;   // Incorrectly oriented fault cells
  PointSet facesSeen;      // Fault faces already considered
  PointSet cellsSeen;      // Fault cells already matched
  Obj<PointSet> newCells  = new PointSet();
  Obj<PointSet> loopCells = new PointSet();

  if (!(*fault)->commRank()) {
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
      const Obj<sieve_type::traits::coneSequence>&     cone   = faultSieve->cone(*c_iter);
      const sieve_type::traits::coneSequence::iterator eBegin = cone->begin();
      const sieve_type::traits::coneSequence::iterator eEnd   = cone->end();

      for(sieve_type::traits::coneSequence::iterator e_iter = eBegin; e_iter != eEnd; ++e_iter) {
        if (facesSeen.find(*e_iter) != facesSeen.end()) continue;
        facesSeen.insert(*e_iter);
        if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
        const Obj<sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
        sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

        // Throw out boundary fault faces
        if (support->size() < 2) continue;
        Mesh::point_type cellA = *s_iter; ++s_iter;
        Mesh::point_type cellB = *s_iter;
        bool flippedA = (flippedCells.find(cellA) != flippedCells.end());
        bool flippedB = (flippedCells.find(cellB) != flippedCells.end());
        bool seenA    = (cellsSeen.find(cellA) != cellsSeen.end());
        bool seenB    = (cellsSeen.find(cellB) != cellsSeen.end());

        if (!seenA) newCells->insert(cellA);
        if (!seenB) newCells->insert(cellB);
        if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
        // In 1D, just check that vertices match
        if (fDim == 1) {
          const Obj<sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
          sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
          const Obj<sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
          sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
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
          ALE::MinimalArrow<sieve_type::point_type,sieve_type::point_type> arrowA(*e_iter, cellA);
          const int oA = orientation->restrictPoint(arrowA)[0];
          ALE::MinimalArrow<sieve_type::point_type,sieve_type::point_type> arrowB(*e_iter, cellB);
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
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = faultSieve->cone(*f_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin(); v_iter != cone->end(); ++v_iter) {
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
        ALE::MinimalArrow<sieve_type::point_type,sieve_type::point_type> arrow(*e_iter, *f_iter);
        int o = orientation->restrictPoint(arrow)[0];

        if (debug) std::cout << "    Reversing orientation of " << *e_iter <<"-->"<<*f_iter << " from " << o << " to " << -(o+1) << std::endl;
        o = -(o+1);
        orientation->updatePoint(arrow, &o);
      }
    }
  }
  flippedCells.clear();
  for(Mesh::label_sequence::iterator e_iter = fFaces->begin(); e_iter != fFaces->end(); ++e_iter) {
    if (debug) std::cout << "  Checking orientation of fault face " << *e_iter << std::endl;
    // for each face get the support (2 fault cells)
    const Obj<sieve_type::traits::supportSequence>& support = faultSieve->support(*e_iter);
    sieve_type::traits::supportSequence::iterator   s_iter  = support->begin();

    // Throw out boundary fault faces
    if (support->size() > 1) {
      Mesh::point_type cellA = *s_iter; ++s_iter;
      Mesh::point_type cellB = *s_iter;

      if (debug) std::cout << "    neighboring cells " << cellA << " and " << cellB << std::endl;
      // In 1D, just check that vertices match
      if (fDim == 1) {
        const Obj<sieve_type::traits::coneSequence>& coneA = faultSieve->cone(cellA);
        sieve_type::traits::coneSequence::iterator   iterA = coneA->begin();
        const Obj<sieve_type::traits::coneSequence>& coneB = faultSieve->cone(cellB);
        sieve_type::traits::coneSequence::iterator   iterB = coneB->begin();
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
        ALE::MinimalArrow<sieve_type::point_type,sieve_type::point_type> arrowA(*e_iter, cellA);
        const int oA = orientation->restrictPoint(arrowA)[0];
        ALE::MinimalArrow<sieve_type::point_type,sieve_type::point_type> arrowB(*e_iter, cellB);
        const int oB = orientation->restrictPoint(arrowB)[0];

        if (oA == oB) {
          std::cout << "Invalid orientation in fault mesh" << std::endl;
          std::cout << "  fault face: " << *e_iter << "  cellA: " << cellA << "  cellB: " << cellB << std::endl;
          throw ALE::Exception("Invalid orientation in fault mesh");
        }
      }
    }
  }
  if (debug) (*fault)->view("Oriented Fault mesh");

  // Add new shadow vertices and possibly Lagrange multipler vertices
  const ALE::Obj<Mesh::label_sequence>& fVertices = (*fault)->depthStratum(0);
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<std::set<std::string> >& groupNames = mesh->getIntSections();
  Mesh::point_type newPoint = sieve->base()->size() + sieve->cap()->size();
  const int        numFaultVertices = fVertices->size();
  std::map<int,int> vertexRenumber;

  for(Mesh::label_sequence::iterator v_iter = fVertices->begin();
      v_iter != fVertices->end();
      ++v_iter, ++newPoint) {
    vertexRenumber[*v_iter] = newPoint;
    if (debug) 
      std::cout << "Duplicating " << *v_iter << " to "
		<< vertexRenumber[*v_iter] << std::endl;

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
      if (group->hasPoint(*v_iter))
        group->addPoint(newPoint, 1);
    } // for
  } // for
  for(std::set<std::string>::const_iterator name = groupNames->begin();
      name != groupNames->end(); ++name) {
    mesh->reallocate(mesh->getIntSection(*name));
  } // for
  if (constraintCell) newPoint += numFaultVertices;

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<Mesh::label_sequence>& faces = (*fault)->heightStratum(1);
  const ALE::Obj<Mesh::label_type>& material = mesh->getLabel("material-id");
  const int firstCohesiveCell = newPoint;
  PointSet replaceCells;
  PointSet noReplaceCells;
  PointSet replaceVertices;

  for(Mesh::label_sequence::iterator f_iter = faces->begin();
      f_iter != faces->end(); ++f_iter, ++newPoint) {
    if (debug) std::cout << "Considering fault face " << *f_iter << std::endl;
    const ALE::Obj<sieve_type::traits::supportSequence>& cells =
      faultSieve->support(*f_iter);
    Mesh::point_type cell = *cells->begin();
    Mesh::point_type otherCell;

    if (debug) std::cout << "  Checking orientation against cell " << cell << std::endl;
    selection::getOrientedFace(mesh, cell, &vertexRenumber, numCorners, indices, &origVertices, &faceVertices);

    const ALE::Obj<sieve_type::coneArray>& faceCone = sieveAlg::nCone(*fault, *f_iter, faultDepth);
    const sieve_type::coneArray::iterator  fBegin   = faceCone->begin();
    const sieve_type::coneArray::iterator  fEnd     = faceCone->end();
    bool found = true;

    if (numFaultCorners == 0) {
      found = false;
    } else if (numFaultCorners == 2) {
      if (faceVertices[0] != *fBegin) found = false;
    } else {
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVertices[v] != *fBegin)) ++v;
      for(sieve_type::coneArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter, ++v) {
        if (debug) std::cout << "    Checking " << *v_iter << " against " << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != *v_iter) {
          found = false;
          break;
        }
      }
    }
    if (found) {
      if (debug) std::cout << "  Choosing other cell" << std::endl;
      otherCell = cell;
      cell = *(++cells->begin());
    } else {
      otherCell = *(++cells->begin());
      if (debug) std::cout << "  Verifing reverse orientation" << std::endl;
      found = true;
      int v = 0;
      if (numFaultCorners > 0) {
        // Locate first vertex
        while((v < numFaultCorners) && (faceVertices[v] != *faceCone->rbegin())) ++v;
        for(sieve_type::coneArray::reverse_iterator v_iter = faceCone->rbegin(); v_iter != faceCone->rend(); ++v_iter, ++v) {
          if (debug) std::cout << "    Checking " << *v_iter << " against " << faceVertices[v%numFaultCorners] << std::endl;
          if (faceVertices[v%numFaultCorners] != *v_iter) {
            found = false;
            break;
          }
        }
      }
      assert(found);
    }
    noReplaceCells.insert(otherCell);
    replaceCells.insert(cell);
    replaceVertices.insert(fBegin, fEnd);
    // Adding cohesive cell (not interpolated)
    int color = 0;

	if (debug)
	  std::cout << "  Creating cohesive cell " << newPoint << std::endl;
    //for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
    for(sieve_type::coneArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    vertex " << *v_iter << std::endl;
      sieve->addArrow(*v_iter, newPoint, color++);
    }
    //for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
    for(sieve_type::coneArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    shadow vertex " << vertexRenumber[*v_iter] << std::endl;
      sieve->addArrow(vertexRenumber[*v_iter], newPoint, color++);
    }
    if (constraintCell) {
      for(sieve_type::coneArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
        if (debug)
          std::cout << "    Lagrange vertex " << vertexRenumber[*v_iter]+numFaultVertices << std::endl;
        sieve->addArrow(vertexRenumber[*v_iter]+numFaultVertices, newPoint, color++);
      }
    }
    mesh->setValue(material, newPoint, materialId);
  } // for
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
    replacedCells->setFiberDimension(mesh->heightStratum(0), 1);
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
  for(PointSet::const_iterator c_iter = replaceCells.begin(); c_iter != replaceCells.end(); ++c_iter) {
    _replaceCell(sieve, *c_iter, &vertexRenumber, debug);
  }
  if (!(*fault)->commRank()) delete [] indices;
  mesh->stratify();
  const std::string labelName("censored depth");

  if (!mesh->hasLabel(labelName)) {
    const ALE::Obj<Mesh::label_type>& label          = mesh->createLabel(labelName);
    const ALE::Obj<PointSet>          modifiedPoints = new PointSet();

    _computeCensoredDepth(mesh, label, mesh->getSieve(), mesh->getSieve()->roots(), firstCohesiveCell-(constraintCell?numFaultVertices:0), modifiedPoints);
  } else {
    // Insert new shadow vertices into existing label
    const ALE::Obj<Mesh::label_type>& label          = mesh->getLabel(labelName);

    for(std::map<int,int>::const_iterator v_iter = vertexRenumber.begin(); v_iter != vertexRenumber.end(); ++v_iter) {
      mesh->setValue(label, v_iter->second, 0);
    }
  }
  if (debug) mesh->view("Mesh with Cohesive Elements");

  // Fix coordinates
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& fVertices2 = (*fault)->depthStratum(0);

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
} // createCohesiveCells

// ----------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::CohesiveTopology::createParallel(
		ALE::Obj<Mesh>* fault,
		std::map<Mesh::point_type, Mesh::point_type>* cohesiveToFault,
		const ALE::Obj<Mesh>& mesh,
		const int materialId,
		const bool constraintCell)
{
  assert(0 != fault);
  assert(0 != cohesiveToFault);

  *fault = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  cohesiveToFault->clear();

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<sieve_type> faultSieve = new sieve_type(sieve->comm(), sieve->debug());
  const ALE::Obj<Mesh::label_sequence>& cohesiveCells = mesh->getLabelStratum("material-id", materialId);
  const Mesh::label_sequence::iterator cBegin = cohesiveCells->begin();
  const Mesh::label_sequence::iterator cEnd = cohesiveCells->end();
  const int sieveEnd = sieve->base()->size() + sieve->cap()->size();
  const int numFaces = cohesiveCells->size();
  int globalSieveEnd = 0;
  int globalFaceOffset = 0;

  MPI_Allreduce((void *) &sieveEnd, (void *) &globalSieveEnd, 1, MPI_INT, MPI_SUM, sieve->comm());
  MPI_Scan((void *) &numFaces, (void *) &globalFaceOffset, 1, MPI_INT, MPI_SUM, sieve->comm());
  int face = globalSieveEnd + globalFaceOffset - numFaces;
  for(Mesh::label_sequence::iterator c_iter = cBegin; c_iter != cEnd; ++c_iter) {
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = sieve->cone(*c_iter);
    const int coneSize = cone->size();
    int       color    = 0;

    if (!constraintCell) {
      const int faceSize = coneSize / 2;
      assert(0 == coneSize % faceSize);

      // Use first vertices (negative side of the fault) for fault mesh
      sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
      for(int i=0; i < faceSize; ++i, ++v_iter) {
        faultSieve->addArrow(*v_iter, face, color++);
      }
    } else {
      const int faceSize = coneSize / 3;
      assert(0 == coneSize % faceSize);

      // Use last vertices (contraints) for fault mesh
      sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
      for(int i=0; i < 2*faceSize; ++i) ++v_iter;
      for(int i=0; i < faceSize; ++i, ++v_iter) {
        faultSieve->addArrow(*v_iter, face, color++);
      }
    } // if/else
    (*cohesiveToFault)[*c_iter] = face;
    ++face;
  } // for
  (*fault)->setSieve(faultSieve);
  (*fault)->stratify();

  const ALE::Obj<Mesh::label_sequence>& faultCells = (*fault)->heightStratum(0);
  assert(!faultCells.isNull());
  for(Mesh::label_sequence::iterator c_iter = cBegin, f_iter=faultCells->begin(); c_iter != cEnd; ++c_iter, ++f_iter) {
    (*cohesiveToFault)[*c_iter] = *f_iter;
  }
    
#if 1
  (*fault)->setRealSection("coordinates", mesh->getRealSection("coordinates"));
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
  const ALE::Obj<sieve_type::traits::supportSequence>& neighbors = sieve->support(vertex);
  const sieve_type::traits::supportSequence::iterator  begin     = neighbors->begin();
  const sieve_type::traits::supportSequence::iterator  end       = neighbors->end();
  bool                                                 modified  = true;
  int                                                  classifyTotal = neighbors->size();
  int                                                  classifySize  = 0;
  PointSet                                             vReplaceCells;
  PointSet                                             vNoReplaceCells;


  if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
  for(sieve_type::traits::supportSequence::iterator n_iter = begin; n_iter != end; ++n_iter) {
    if (replaceCells.find(*n_iter)   != replaceCells.end())   vReplaceCells.insert(*n_iter);
    if (noReplaceCells.find(*n_iter) != noReplaceCells.end()) vNoReplaceCells.insert(*n_iter);
    if (*n_iter >= firstCohesiveCell) classifyTotal--;
  }
  classifySize = vReplaceCells.size() + vNoReplaceCells.size();

  while(modified && (classifySize < classifyTotal)) {
    modified = false;
    for(sieve_type::traits::supportSequence::iterator n_iter = begin; n_iter != end; ++n_iter) {
      bool classified = false;

      if (debug) {std::cout << "Checking neighbor " << *n_iter << std::endl;}
      if (vReplaceCells.find(*n_iter)   != vReplaceCells.end()) {
        if (debug) {std::cout << "  already in replaceCells" << std::endl;}
        continue;
      }
      if (vNoReplaceCells.find(*n_iter) != vNoReplaceCells.end()) {
        if (debug) {std::cout << "  already in noReplaceCells" << std::endl;}
        continue;
      }
      if (*n_iter >= firstCohesiveCell) {
        if (debug) {std::cout << "  already a cohesive cell" << std::endl;}
        continue;
      }
      // If neighbor shares a face with anyone in replaceCells, then add
      for(PointSet::const_iterator c_iter = vReplaceCells.begin(); c_iter != vReplaceCells.end(); ++c_iter) {
        const ALE::Obj<sieve_type::coneSet>& preFace = sieve->nMeet(*c_iter, *n_iter, depth);

        if (preFace->size() == faceSize) {
          if (debug) {std::cout << "    Scheduling " << *n_iter << " for replacement" << std::endl;}
          vReplaceCells.insert(*n_iter);
          modified   = true;
          classified = true;
          break;
        }
      }
      if (classified) continue;
      // It is unclear whether taking out the noReplace cells will speed this up
      for(PointSet::const_iterator c_iter = vNoReplaceCells.begin(); c_iter != vNoReplaceCells.end(); ++c_iter) {
        const ALE::Obj<sieve_type::coneSet>& preFace = sieve->nMeet(*c_iter, *n_iter, depth);

        if (preFace->size() == faceSize) {
          if (debug) {std::cout << "    Scheduling " << *n_iter << " for no replacement" << std::endl;}
          vNoReplaceCells.insert(*n_iter);
          modified   = true;
          classified = true;
          break;
        }
      }
    }
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
pylith::faults::CohesiveTopology::_replaceCell(const Obj<sieve_type>& sieve,
                                               const Mesh::point_type cell,
                                               std::map<int,int> *vertexRenumber,
                                               const int debug)
{
  bool       replace = false;
  PointArray newVertices;

  const ALE::Obj<sieve_type::traits::coneSequence>& cCone = sieve->cone(cell);

  for(sieve_type::traits::coneSequence::iterator v_iter = cCone->begin();
      v_iter != cCone->end(); ++v_iter) {
    if (vertexRenumber->find(*v_iter) != vertexRenumber->end()) {
      if (debug) std::cout << "    vertex " << (*vertexRenumber)[*v_iter] << std::endl;
      newVertices.insert(newVertices.end(), (*vertexRenumber)[*v_iter]);
      replace = true;
    } else {
      if (debug) std::cout << "    vertex " << *v_iter << std::endl;
      newVertices.insert(newVertices.end(), *v_iter);
    } // if/else
  } // for
  if (replace) {
    if (debug) std::cout << "  Replacing cell " << cell << std::endl;
    sieve->clearCone(cell);
    int color = 0;
    for(PointArray::const_iterator v_iter = newVertices.begin(); v_iter != newVertices.end(); ++v_iter) {
      sieve->addArrow(*v_iter, cell, color++);
    } // for
  }
}

// ----------------------------------------------------------------------
template<class InputPoints>
void
pylith::faults::CohesiveTopology::_computeCensoredDepth(const ALE::Obj<Mesh>& mesh,
                                                        const ALE::Obj<Mesh::label_type>& depth,
                                                        const ALE::Obj<Mesh::sieve_type>& sieve,
                                                        const ALE::Obj<InputPoints>& points,
                                                        const Mesh::point_type& firstCohesiveCell,
                                                        const ALE::Obj<std::set<Mesh::point_type> >& modifiedPoints)
{
  modifiedPoints->clear();

  for(typename InputPoints::iterator p_iter = points->begin(); p_iter != points->end(); ++p_iter) {
    if (*p_iter >= firstCohesiveCell) continue;
    // Compute the max depth of the points in the cone of p, and add 1
    int d0 = mesh->getValue(depth, *p_iter, -1);
    int d1 = mesh->getMaxValue(depth, sieve->cone(*p_iter), -1) + 1;

    if(d1 != d0) {
      mesh->setValue(depth, *p_iter, d1);
      modifiedPoints->insert(*p_iter);
    }
  }
  // FIX: We would like to avoid the copy here with support()
  if(modifiedPoints->size() > 0) {
    _computeCensoredDepth(mesh, depth, sieve, sieve->support(modifiedPoints), firstCohesiveCell, modifiedPoints);
  }
}


// End of file
