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

  typedef ALE::SieveAlg<Mesh> sieveAlg;
  typedef std::set<Mesh::point_type> PointSet;

  const int_section_type::chart_type& chart = groupField->getChart();
  PointSet faultVertices; // Vertices on fault
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  *fault = new Mesh(mesh->comm(), mesh->getDimension()-1, mesh->debug());
  const ALE::Obj<sieve_type> faultSieve = new sieve_type(sieve->comm(), 
                                                         sieve->debug());
  const int  numCells   = mesh->heightStratum(0)->size();
  int        numCorners = 0;    // The number of vertices in a mesh cell
  int        faceSize   = 0;    // The number of vertices in a mesh face
  int       *indices    = NULL; // The indices of a face vertex set in a cell
  int        oppositeVertex;    // For simplices, the vertex opposite a given face
  PointArray origVertices;
  PointArray faceVertices;
  PointArray neighborVertices;

  if (!(*fault)->commRank()) {
    numCorners = sieve->nCone(*mesh->heightStratum(0)->begin(), mesh->depth())->size();
    faceSize   = _numFaceVertices(*mesh->heightStratum(0)->begin(), mesh);
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
  for(PointSet::const_iterator fv_iter = fvBegin; fv_iter != fvEnd; ++fv_iter) {
    const ALE::Obj<sieveAlg::supportArray>& cells =
      sieveAlg::nSupport(mesh, *fv_iter, mesh->depth());
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
        const ALE::Obj<sieve_type::supportSet> preFace = faultSieve->nJoin1(face);

        if (preFace->size() > 1) {
          throw ALE::Exception("Invalid fault sieve: Multiple faces from "
                               "vertex set");
        } else if (preFace->size() == 1) {
          // Add the other cell neighbor for this face
          faultSieve->addArrow(*preFace->begin(), *c_iter);
        } else if (preFace->size() == 0) {
          if (debug) std::cout << "  Orienting face " << f << std::endl;
          _getOrientedFace(mesh, *c_iter, face, numCorners, indices, &origVertices, &faceVertices);
          if (debug) std::cout << "  Adding face " << f << std::endl;
          int color = 0;
          for(PointArray::const_iterator f_iter = faceVertices.begin();
              f_iter != faceVertices.end(); ++f_iter) {
            if (debug) std::cout << "    vertex " << *f_iter << std::endl;
            faultSieve->addArrow(*f_iter, f, color++);
          } // for
          faultSieve->addArrow(f, *c_iter);
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
  // Orient the fault sieve
  const ALE::Obj<Mesh::label_sequence>& fFaces = (*fault)->heightStratum(1);
  int faultDepth      = (*fault)->depth()-1; // Depth of fault cells
  int numFaultCorners = 0; // The number of vertices in a fault cell
  int faultFaceSize   = 0; // The number of vertices in a face between fault cells
  PointArray flippedFaces; // Incorrectly oriented fault cells

  if (!(*fault)->commRank()) {
    numFaultCorners = faultSieve->nCone(*fFaces->begin(), faultDepth)->size();
    if (debug) std::cout << "  Fault corners " << numFaultCorners << std::endl;
    assert(numFaultCorners == faceSize);
    faultFaceSize = _numFaceVertices(*fFaces->begin(), (*fault), faultDepth);
  }
  if (debug) std::cout << "  Fault face size " << faultFaceSize << std::endl;
  for(Mesh::label_sequence::iterator e_iter = fFaces->begin(); e_iter != fFaces->end(); ++e_iter) {
    if (debug) std::cout << "  Checking fault face " << *e_iter << std::endl;
    const Obj<sieve_type::traits::coneSequence>& vertices = faultSieve->cone(*e_iter);
    sieve_type::traits::coneSequence::iterator   vEnd     = vertices->end();
    PointSet facesSeen;

    for(sieve_type::traits::coneSequence::iterator v_iter = vertices->begin(); v_iter != vEnd; ++v_iter) {
      const Obj<sieve_type::traits::supportSequence>& neighbors  = faultSieve->support(*v_iter);
      const sieve_type::traits::supportSequence::iterator nBegin = neighbors->begin();
      const sieve_type::traits::supportSequence::iterator nEnd   = neighbors->end();

      for(sieve_type::traits::supportSequence::iterator n_iter = nBegin; n_iter != nEnd; ++n_iter) {
        if (facesSeen.find(*n_iter) != facesSeen.end()) continue;
        facesSeen.insert(*n_iter);
        if (debug) std::cout << "  Checking fault neighbor " << *n_iter << std::endl;
        if (*e_iter >= *n_iter) continue;
        const ALE::Obj<sieve_type::coneSet>& meet = faultSieve->nMeet(*e_iter, *n_iter, faultDepth);

        if (debug) {
          for(sieve_type::coneSet::iterator c_iter = meet->begin(); c_iter != meet->end(); ++c_iter)
            std::cout << "    meet " << *c_iter << std::endl;
        }
        if ((int) meet->size() == faultFaceSize) {
          if (debug) std::cout << "    Found neighboring fault face " << *n_iter << std::endl;
          bool eOrient = _getOrientedFace(*fault, *e_iter, meet, numFaultCorners, indices, &origVertices, &faceVertices);
          bool nOrient = _getOrientedFace(*fault, *n_iter, meet, numFaultCorners, indices, &origVertices, &neighborVertices);

          if (faultFaceSize > 1) {
            if (debug) {
              for(PointArray::iterator v_iter = faceVertices.begin(); v_iter != faceVertices.end(); ++v_iter) {
                std::cout << "  face vertex " << *v_iter << std::endl;
              }
              for(PointArray::iterator v_iter = neighborVertices.begin(); v_iter != neighborVertices.end(); ++v_iter) {
                std::cout << "  neighbor vertex " << *v_iter << std::endl;
              }
            }
            // Here we use the fact that fault faces are only 1D
            if (*faceVertices.begin() == *neighborVertices.begin()) {
              if (debug) std::cout << "  Scheduling fault face " << *n_iter << " to be flipped" << std::endl;
              flippedFaces.push_back(*n_iter);
            }
          } else {
            // For 0D, we use the orientation returned (not sure if we have to do this)
            if (nOrient == eOrient) {
              if (debug) std::cout << "  Scheduling fault face " << *n_iter << " to be flipped" << std::endl;
              flippedFaces.push_back(*n_iter);
            }
          }
        }
      }
    }
  }
  for(PointArray::const_iterator f_iter = flippedFaces.begin(); f_iter != flippedFaces.end(); ++f_iter) {
    if (debug) std::cout << "  Reversing fault face " << *f_iter << std::endl;
    faceVertices.clear();
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = faultSieve->cone(*f_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
        v_iter != cone->end(); ++v_iter) {
      faceVertices.insert(faceVertices.begin(), *v_iter);
    }
    faultSieve->clearCone(*f_iter);
    int color = 0;
    for(PointArray::const_iterator v_iter = faceVertices.begin();
        v_iter != faceVertices.end(); ++v_iter) {
      faultSieve->addArrow(*v_iter, *f_iter, color++);
    } // for
  }
  flippedFaces.clear();
  if (debug) (*fault)->view("Oriented Fault mesh");

  // Add new shadow vertices and possibly Lagrange multipler vertices
  const ALE::Obj<Mesh::label_sequence>& fVertices = (*fault)->depthStratum(0);
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<std::set<std::string> >& groupNames = mesh->getIntSections();
  Mesh::point_type newPoint = sieve->base()->size() + sieve->cap()->size();
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
      groupField->addPoint(newPoint+1, 1);

    // Add shadow vertices to other groups, don't add constraint
    // vertices (if they exist) because we don't want BC, etc to act
    // on constraint vertices
    for(std::set<std::string>::const_iterator name = groupNames->begin();
       name != groupNames->end(); ++name) {
      const ALE::Obj<int_section_type>& group = mesh->getIntSection(*name);
      if (group->hasPoint(*v_iter))
        group->addPoint(newPoint, 1);
    } // for
    if (constraintCell) newPoint++;
  } // for
  for(std::set<std::string>::const_iterator name = groupNames->begin();
      name != groupNames->end(); ++name) {
    mesh->reallocate(mesh->getIntSection(*name));
  } // for

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<Mesh::label_sequence>& faces = (*fault)->depthStratum(1);
  const ALE::Obj<Mesh::label_type>& material = mesh->getLabel("material-id");
  int firstCohesiveCell = newPoint;
  PointArray newVertices;

  for(Mesh::label_sequence::iterator f_iter = faces->begin();
      f_iter != faces->end(); ++f_iter, ++newPoint) {
    if (debug) std::cout << "Considering fault face " << *f_iter << std::endl;
    const ALE::Obj<sieve_type::traits::supportSequence>& cells =
      faultSieve->support(*f_iter);
    Mesh::point_type cell = *cells->begin();

    if (debug) std::cout << "  Checking orientation against cell " << cell << std::endl;
    _getOrientedFace(mesh, cell, &vertexRenumber, numCorners, indices, &origVertices, &faceVertices);

    const ALE::Obj<sieve_type::traits::coneSequence>& faceCone = faultSieve->cone(*f_iter);
    bool found = true;

    if (numFaultCorners == 2) {
      if (faceVertices[0] != *faceCone->begin()) found = false;
    } else {
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVertices[v] != *faceCone->begin())) ++v;
      for(sieve_type::traits::coneSequence::iterator v_iter = faceCone->begin(); v_iter != faceCone->end(); ++v_iter, ++v) {
        if (debug) std::cout << "    Checking " << *v_iter << " against " << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != *v_iter) {
          found = false;
          break;
        }
      }
    }
    if (found) {
      if (debug) std::cout << "  Choosing other cell" << std::endl;
      cell = *(++cells->begin());
    } else {
      if (debug) std::cout << "  Verifing reverse orientation" << std::endl;
      found = true;
      int v = 0;
      // Locate first vertex
      while((v < numFaultCorners) && (faceVertices[v] != *faceCone->rbegin())) ++v;
      for(sieve_type::traits::coneSequence::reverse_iterator v_iter = faceCone->rbegin(); v_iter != faceCone->rend(); ++v_iter, ++v) {
        if (debug) std::cout << "    Checking " << *v_iter << " against " << faceVertices[v%numFaultCorners] << std::endl;
        if (faceVertices[v%numFaultCorners] != *v_iter) {
          found = false;
          break;
        }
      }
      assert(found);
    }
    if (debug) std::cout << "  Replacing cell " << cell << std::endl;
    const ALE::Obj<sieve_type::traits::coneSequence>& cCone = sieve->cone(cell);

    newVertices.clear();
    for(sieve_type::traits::coneSequence::iterator v_iter = cCone->begin();
        v_iter != cCone->end(); ++v_iter) {
      if (vertexRenumber.find(*v_iter) != vertexRenumber.end()) {
        if (debug)
          std::cout << "    vertex " << vertexRenumber[*v_iter] << std::endl;
        newVertices.insert(newVertices.end(), vertexRenumber[*v_iter]);
      } else {
        if (debug) std::cout << "    vertex " << *v_iter << std::endl;
        newVertices.insert(newVertices.end(), *v_iter);
      } // if/else
    } // for
    sieve->clearCone(cell);
    int color = 0;
    for(PointArray::const_iterator v_iter = newVertices.begin();
        v_iter != newVertices.end(); ++v_iter) {
      sieve->addArrow(*v_iter, cell, color++);
    } // for
    // Adding cohesive cell (not interpolated)
    const ALE::Obj<sieve_type::traits::coneSequence>& fCone  = faultSieve->cone(*f_iter);
    const sieve_type::traits::coneSequence::iterator  fBegin = fCone->begin();
    const sieve_type::traits::coneSequence::iterator  fEnd   = fCone->end();
    color = 0;

	if (debug)
	  std::cout << "  Creating cohesive cell " << newPoint << std::endl;
    //for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
    for(sieve_type::traits::coneSequence::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    vertex " << *v_iter << std::endl;
      sieve->addArrow(*v_iter, newPoint, color++);
    }
    //for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
    for(sieve_type::traits::coneSequence::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    shadow vertex " << vertexRenumber[*v_iter] << std::endl;
      sieve->addArrow(vertexRenumber[*v_iter], newPoint, color++);
    }
    if (constraintCell) {
      //for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      for(sieve_type::traits::coneSequence::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
        if (debug)
          std::cout << "    Lagrange vertex " << vertexRenumber[*v_iter]+1 << std::endl;
        sieve->addArrow(vertexRenumber[*v_iter]+1, newPoint, color++);
      }
    }
    mesh->setValue(material, newPoint, materialId);
  } // for
  if (!(*fault)->commRank()) delete [] indices;
  mesh->stratify();
  const ALE::Obj<Mesh::label_type>& label          = mesh->createLabel(std::string("censored depth"));
  const ALE::Obj<PointSet>          modifiedPoints = new PointSet();
  _computeCensoredDepth(mesh, label, mesh->getSieve(), mesh->getSieve()->roots(), firstCohesiveCell, modifiedPoints);
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
      coordinates->addPoint(vertexRenumber[*v_iter]+1,
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
      coordinates->updatePoint(vertexRenumber[*v_iter]+1,
			     coordinates->restrictPoint(*v_iter));
    }
  }
} // createCohesiveCells

// ----------------------------------------------------------------------
unsigned int
pylith::faults::CohesiveTopology::_numFaceVertices(const Mesh::point_type& cell,
                                                   const ALE::Obj<Mesh>& mesh,
                                                   const int depth)
{ // _numFaceVertices

  /** :QUESTION:
   *
   * If mesh is interpolated is there a simple way to get the number
   * of vertices on the face (3-D), edge (2-D), end (1-D) of a cell?
   */
  const int cellDim = mesh->getDimension();

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  int meshDepth = (depth < 0) ? mesh->depth() : depth;
  unsigned int numCorners = sieve->nCone(cell, meshDepth)->size();

  unsigned int numFaceVertices = 0;
  switch (cellDim)
    { // switch
    case 0 :
      numFaceVertices = 0;
      break;
    case 1 :
      numFaceVertices = 1;
      break;
    case 2:
      switch (numCorners)
	{ // switch
	case 3 : // tri3
	  numFaceVertices = 2; // Edge has 2 vertices
	  break;
	case 4 : // quad4
	  numFaceVertices = 2; // Edge has 2 vertices
	  break;
	default :
	  std::cerr << "numCorners: " << numCorners << std::endl;
	  assert(0);
	} // switch
      break;
    case 3:
      switch (numCorners)
	{ // switch
	case 4 : // tet4
	  numFaceVertices = 3; // Face has 3 vertices
	  break;
	case 8 : // hex8
	  numFaceVertices = 4; // Face has 4 vertices
	  break;
	default :
	  std::cerr << "numCorners: " << numCorners << std::endl;
	  assert(0);
	} // switch
      break;
    default:
      assert(0);
    } // swtich
  return numFaceVertices;
} // _numFaceVertices

// ----------------------------------------------------------------------
// We need this method because we do not use interpolates sieves
//   - Without interpolation, we cannot say what vertex collections are
//     faces, and how they are oriented
//   - Now we read off the list of face vertices IN THE ORDER IN WHICH
//     THEY APPEAR IN THE CELL
//   - This leads to simple algorithms for simplices and quads to check
//     orientation since these sets are always valid faces
//   - This is not true with hexes, so we just sort and check explicit cases
//   - This means for hexes that we have to alter the vertex container as well
bool
pylith::faults::CohesiveTopology::_faceOrientation(const Mesh::point_type& cell,
                                                   const ALE::Obj<Mesh>& mesh,
                                                   const int numCorners,
                                                   const int indices[],
                                                   const int oppositeVertex,
                                                   PointArray *origVertices,
                                                   PointArray *faceVertices)
{ // _faceOrientation
  const int cellDim   = mesh->getDimension();
  bool      posOrient = false;
  int       debug     = mesh->debug();

  // Simplices
  if (cellDim == numCorners-1) {
    posOrient = !(oppositeVertex%2);
  } else if (cellDim == 2) {
    // Quads
    if ((indices[1] > indices[0]) && (indices[1] - indices[0] == 1)) {
      posOrient = true;
    } else if ((indices[0] == 3) && (indices[1] == 0)) {
      posOrient = true;
    } else {
      posOrient = false;
    }
  } else if (cellDim == 3) {
    // Hexes
    //   A hex is two oriented quads with the normal of the first
    //   pointing up at the second.
    //
    //     7---6
    //    /|  /|
    //   4---5 |
    //   | 3-|-2
    //   |/  |/
    //   0---1
    int sortedIndices[4];

    for(int i = 0; i < 4; ++i) sortedIndices[i] = indices[i];
    std::sort(sortedIndices, sortedIndices+4);
    // Case 1: Bottom quad
    if ((sortedIndices[0] == 0) && (sortedIndices[1] == 1) && (sortedIndices[2] == 2) && (sortedIndices[3] == 3)) {
      if (debug) std::cout << "Bottom quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 3) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 2) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 1) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 0) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    // Case 2: Top quad
    if ((sortedIndices[0] == 4) && (sortedIndices[1] == 5) && (sortedIndices[2] == 6) && (sortedIndices[3] == 7)) {
      if (debug) std::cout << "Top quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 5) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 6) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 7) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 4) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    // Case 3: Front quad
    if ((sortedIndices[0] == 0) && (sortedIndices[1] == 1) && (sortedIndices[2] == 4) && (sortedIndices[3] == 5)) {
      if (debug) std::cout << "Front quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 1) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 5) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 4) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 0) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    // Case 4: Back quad
    if ((sortedIndices[0] == 2) && (sortedIndices[1] == 3) && (sortedIndices[2] == 6) && (sortedIndices[3] == 7)) {
      if (debug) std::cout << "Back quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 7) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 6) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 2) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 3) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    // Case 5: Right quad
    if ((sortedIndices[0] == 1) && (sortedIndices[1] == 2) && (sortedIndices[2] == 5) && (sortedIndices[3] == 6)) {
      if (debug) std::cout << "Right quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 2) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 6) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 5) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 1) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    // Case 6: Left quad
    if ((sortedIndices[0] == 0) && (sortedIndices[1] == 3) && (sortedIndices[2] == 4) && (sortedIndices[3] == 7)) {
      if (debug) std::cout << "Left quad" << std::endl;
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 4) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 7) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 3) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
      for(int i = 0; i < 4; ++i) {
        if (indices[i] == 0) {
          faceVertices->push_back((*origVertices)[i]); break;
        }
      }
    }
    return true;
  }
  if (!posOrient) {
    if (debug) std::cout << "  Reversing initial face orientation" << std::endl;
    faceVertices->insert(faceVertices->end(), (*origVertices).rbegin(), (*origVertices).rend());
  } else {
    if (debug) std::cout << "  Keeping initial face orientation" << std::endl;
    faceVertices->insert(faceVertices->end(), (*origVertices).begin(), (*origVertices).end());
  }
  return posOrient;
} // _faceOrientation

// ----------------------------------------------------------------------
// Given a cell and a face, as a set of vertices,
//   return the oriented face, as a set of vertices, in faceVertices
// The orientation is such that the face normal points out of the cell
template<typename FaceType>
bool
pylith::faults::CohesiveTopology::_getOrientedFace(const ALE::Obj<Mesh>& mesh,
                                                   const Mesh::point_type& cell,
                                                   FaceType face,
                                                   const int numCorners,
                                                   int indices[],
                                                   PointArray *origVertices,
                                                   PointArray *faceVertices)
{
  const ALE::Obj<sieve_type::traits::coneSequence>& cone = mesh->getSieve()->cone(cell);
  const int debug = mesh->debug();
  int       v     = 0;
  int       oppositeVertex;

  origVertices->clear();
  faceVertices->clear();
  for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
      v_iter != cone->end(); ++v_iter, ++v) {
    if (face->find(*v_iter) != face->end()) {
      if (debug) std::cout << "    vertex " << *v_iter << std::endl;
      indices[origVertices->size()] = v;
      origVertices->insert(origVertices->end(), *v_iter);
    } else {
      if (debug) std::cout << "    vertex " << *v_iter << std::endl;
      oppositeVertex = v;
    }
  }
  return _faceOrientation(cell, mesh, numCorners, indices, oppositeVertex, origVertices, faceVertices);
}

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
};


// End of file
