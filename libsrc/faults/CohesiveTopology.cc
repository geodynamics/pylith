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

  typedef std::vector<Mesh::point_type> PointArray;
  typedef ALE::SieveAlg<Mesh> sieveAlg;

  // Create set with vertices on fault
  const int_section_type::chart_type& chart = groupField->getChart();
  std::set<Mesh::point_type> faultVertices; // Vertices on fault

  const int numCells = mesh->heightStratum(0)->size();
  for(int_section_type::chart_type::iterator c_iter = chart.begin();
      c_iter != chart.end();
      ++c_iter) {
    assert(!mesh->depth(*c_iter));
    faultVertices.insert(*c_iter);
  } // for

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  *fault = new Mesh(mesh->comm(), mesh->debug());
  const ALE::Obj<sieve_type> faultSieve = new sieve_type(sieve->comm(), 
						    sieve->debug());
  const std::set<Mesh::point_type>::const_iterator fvBegin = 
    faultVertices.begin();
  const std::set<Mesh::point_type>::const_iterator fvEnd = 
    faultVertices.end();

  int f = sieve->base()->size() + sieve->cap()->size();
  int debug = mesh->debug();
  ALE::Obj<PointArray> face = new PointArray();
  std::set<Mesh::point_type> faultCells;
  
  // Create a sieve which captures the fault
  for(std::set<int>::const_iterator fv_iter = fvBegin;
      fv_iter != fvEnd;
      ++fv_iter) {
    const ALE::Obj<sieveAlg::supportArray>& cells =
      sieveAlg::nSupport(mesh, *fv_iter, mesh->depth());
    const sieveAlg::supportArray::iterator cBegin = cells->begin();
    const sieveAlg::supportArray::iterator cEnd   = cells->end();
    
    if (debug)
      std::cout << "Checking fault vertex " << *fv_iter << std::endl;
    for(sieveAlg::supportArray::iterator c_iter = cBegin;
	c_iter != cEnd;
	++c_iter) {
      const unsigned int faceSize = _numFaceVertices(*c_iter, mesh);

      if (debug)
	std::cout << "  Checking cell " << *c_iter << std::endl;
      if (faultCells.find(*c_iter) != faultCells.end())
	continue;
      const ALE::Obj<sieveAlg::coneArray>& cone =
        sieveAlg::nCone(mesh, *c_iter, mesh->height());
      const sieveAlg::coneArray::iterator vBegin = cone->begin();
      const sieveAlg::coneArray::iterator vEnd   = cone->end();
      
      face->clear();
      for(sieveAlg::coneArray::iterator v_iter = vBegin;
	  v_iter != vEnd;
	  ++v_iter) {
	if (faultVertices.find(*v_iter) != fvEnd) {
	  if (debug)
	    std::cout << "    contains fault vertex " << *v_iter << std::endl;
	  face->insert(face->end(), *v_iter);
	} // if
      } // for
      if (face->size() > faceSize)
	throw ALE::Exception("Invalid fault mesh: Too many vertices of an "
			     "element on the fault");
      if (face->size() == faceSize) {
        if (debug)
          std::cout << "  Contains a face on the fault" << std::endl;
        const ALE::Obj<sieve_type::supportSet> preFace = 
          faultSieve->nJoin1(face);
	
        if (preFace->size() > 1)
          throw ALE::Exception("Invalid fault sieve: Multiple faces from "
                               "vertex set");
        else if (preFace->size() == 1)
          faultSieve->addArrow(*preFace->begin(), *c_iter);
        else if (preFace->size() == 0) {
          if (debug)
            std::cout << "  Adding face " << f << std::endl;
          int color = 0;
          for(PointArray::const_iterator f_iter = face->begin();
              f_iter != face->end();
              ++f_iter) {
            if (debug)
              std::cout << "    vertex " << *f_iter << std::endl;
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
  if (debug)
    (*fault)->view("Fault mesh");

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

    for(std::set<std::string>::const_iterator name = groupNames->begin();
       name != groupNames->end(); ++name) {
      const ALE::Obj<int_section_type>& group = mesh->getIntSection(*name);

      if (group->hasPoint(*v_iter)) {
        group->addPoint(newPoint, 1);
        if (constraintCell) {
          group->addPoint(newPoint+1, 1);
        }
      }
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
  PointArray origVertices;
  PointArray newVertices;
  int        oppositeVertex;
  const int  numCorners = sieve->nCone(*mesh->heightStratum(0)->begin(), mesh->depth())->size();
  const int  faceSize   = _numFaceVertices(*mesh->heightStratum(0)->begin(), mesh);
  int       *indices    = new int[faceSize];
  const int  firstCohesiveCell = newPoint;
  
  for(Mesh::label_sequence::iterator f_iter = faces->begin();
      f_iter != faces->end();
      ++f_iter, ++newPoint) {
    if (debug)
      std::cout << "Considering fault face " << *f_iter << std::endl;
    const ALE::Obj<sieve_type::traits::supportSequence>& cells =
      faultSieve->support(*f_iter);
    Mesh::point_type cell = std::max(*cells->begin(), *(++cells->begin()));
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = sieve->cone(cell);
    
    if (debug)
      std::cout << "  Replacing cell " << cell << std::endl;
    origVertices.clear();
    newVertices.clear();
    int v = 0;
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
        v_iter != cone->end();
        ++v_iter, ++v) {
      if (vertexRenumber.find(*v_iter) != vertexRenumber.end()) {
        if (debug)
          std::cout << "    vertex " << vertexRenumber[*v_iter] << std::endl;
        indices[origVertices.size()] = v;
        newVertices.insert(newVertices.end(), vertexRenumber[*v_iter]);
        origVertices.insert(origVertices.end(), *v_iter);
      } else {
        if (debug)
          std::cout << "    vertex " << *v_iter << std::endl;
        newVertices.insert(newVertices.end(), *v_iter);
        oppositeVertex = v;
      } // if/else
    } // for
    sieve->clearCone(cell);
    int color = 0;
    for(PointArray::const_iterator v_iter = newVertices.begin();
        v_iter != newVertices.end();
        ++v_iter) {
      sieve->addArrow(*v_iter, cell, color++);
    } // for
    // Adding cohesive cell (not interpolated)
    //const ALE::Obj<sieve_type::traits::coneSequence>& fCone  = faultSieve->cone(*f_iter);
    //const sieve_type::traits::coneSequence::iterator  fBegin = fCone->begin();
    //const sieve_type::traits::coneSequence::iterator  fEnd   = fCone->end();
    PointArray faceVertices;
	if (debug) {
      int v = 0;

	  std::cout << "  Original Vertices: " << std::endl << "    ";
      for(PointArray::iterator v_iter = origVertices.begin(); v_iter != origVertices.end(); ++v_iter, ++v) {
        std::cout << " " << *v_iter << "(" << indices[v] << ")";
      }
	  std::cout << std::endl << "  Opposite Vertex: " << oppositeVertex << std::endl;
    }
    if (_faceOrientation(cell, mesh, numCorners, indices, oppositeVertex)) {
      if (debug)
        std::cout << "  Reversing initial face orientation" << std::endl;
      faceVertices.insert(faceVertices.end(), origVertices.rbegin(), origVertices.rend());
    } else {
      if (debug)
        std::cout << "  Keeping initial face orientation" << std::endl;
      faceVertices.insert(faceVertices.end(), origVertices.begin(), origVertices.end());
    }
    const PointArray::iterator fBegin = faceVertices.begin();
    const PointArray::iterator fEnd   = faceVertices.end();
    color = 0;

	if (debug)
	  std::cout << "  Creating cohesive cell " << newPoint << std::endl;
    for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    vertex " << *v_iter << std::endl;
      sieve->addArrow(*v_iter, newPoint, color++);
    }
    for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
      if (debug)
        std::cout << "    shadow vertex " << vertexRenumber[*v_iter] << std::endl;
      sieve->addArrow(vertexRenumber[*v_iter], newPoint, color++);
    }
    if (constraintCell) {
      for(PointArray::iterator v_iter = fBegin; v_iter != fEnd; ++v_iter) {
        if (debug)
          std::cout << "    Lagrange vertex " << vertexRenumber[*v_iter]+1 << std::endl;
        sieve->addArrow(vertexRenumber[*v_iter]+1, newPoint, color++);
      }
    }
    mesh->setValue(material, newPoint, materialId);
  } // for
  delete [] indices;
  mesh->stratify();
  const ALE::Obj<Mesh::label_type>&           label          = mesh->createLabel(std::string("censored depth"));
  const ALE::Obj<std::set<Mesh::point_type> > modifiedPoints = new std::set<Mesh::point_type>();
  _computeCensoredDepth(mesh, label, mesh->getSieve(), mesh->getSieve()->roots(), firstCohesiveCell, modifiedPoints);
  if (debug)
    mesh->view("Mesh with Cohesive Elements");

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
						   const ALE::Obj<Mesh>& mesh)
{ // _numFaceVertices

  /** :QUESTION:
   *
   * If mesh is interpolated is there a simple way to get the number
   * of vertices on the face (3-D), edge (2-D), end (1-D) of a cell?
   */
  const int cellDim = mesh->getDimension();

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  unsigned int numCorners = sieve->nCone(cell, mesh->depth())->size();

  unsigned int numFaceVertices = 0;
  switch (cellDim)
    { // switch
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
bool
pylith::faults::CohesiveTopology::_faceOrientation(const Mesh::point_type& cell,
                                                   const ALE::Obj<Mesh>& mesh,
                                                   const int numCorners,
                                                   const int indices[],
                                                   const int oppositeVertex)
{ // _faceOrientation
  const int cellDim = mesh->getDimension();

  // Simplices
  if (cellDim == numCorners-1) {
    return !(oppositeVertex%2);
  } else if (cellDim == 2) {
    // Quads
    if ((indices[1] > indices[0]) && (indices[1] - indices[0] == 1)) {
      return true;
    }
    return false;
  } else if (cellDim == 3) {
    // Hexes
    //   I think we might have to enumerate all of these, ugh
  }
  return true;
} // _faceOrientation

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
