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
pylith::faults::CohesiveTopology::create(const ALE::Obj<Mesh>& mesh,
					 const ALE::Obj<Mesh::int_section_type>& groupField,
					 const int materialId)
{ // create
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
  const ALE::Obj<Mesh> fault = new Mesh(mesh->comm(), mesh->debug());
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
  fault->setSieve(faultSieve);
  fault->stratify();
  faultCells.clear();
  if (debug)
    fault->view("Fault mesh");

  // Add new shadow vertices
  const ALE::Obj<Mesh::label_sequence>& fVertices = fault->depthStratum(0);
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  Mesh::point_type newPoint = sieve->base()->size() + sieve->cap()->size();
  std::map<int,int> vertexRenumber;
  
  for(Mesh::label_sequence::iterator v_iter = fVertices->begin();
      v_iter != fVertices->end();
      ++v_iter, ++newPoint) {
    if (debug) 
      std::cout << "Duplicating " << *v_iter << " to "
		<< vertexRenumber[*v_iter] << std::endl;
    vertexRenumber[*v_iter] = newPoint;
    groupField->setFiberDimension(newPoint, 1);
  } // for

  // Split the mesh along the fault sieve and create cohesive elements
  const ALE::Obj<Mesh::label_sequence>& faces = fault->depthStratum(1);
  const ALE::Obj<Mesh::label_type>& material = mesh->getLabel("material-id");
  PointArray newVertices;
  
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
    newVertices.clear();
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
	v_iter != cone->end();
	++v_iter) {
      if (vertexRenumber.find(*v_iter) != vertexRenumber.end()) {
	if (debug)
	  std::cout << "    vertex " << vertexRenumber[*v_iter] << std::endl;
	newVertices.insert(newVertices.end(), vertexRenumber[*v_iter]);
      } else {
	if (debug)
	  std::cout << "    vertex " << *v_iter << std::endl;
	newVertices.insert(newVertices.end(), *v_iter);
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
    const ALE::Obj<sieve_type::traits::coneSequence>& fCone  = faultSieve->cone(*f_iter);
    const sieve_type::traits::coneSequence::iterator  fBegin = fCone->begin();
    const sieve_type::traits::coneSequence::iterator  fEnd   = fCone->end();
    color = 0;

	if (debug)
	  std::cout << "  Creating cohesive cell " << newPoint << std::endl;
    for(sieve_type::traits::coneSequence::iterator v_iter = fBegin; v_iter != fEnd;
        ++v_iter) {
      if (debug)
        std::cout << "    vertex " << *v_iter << std::endl;
      sieve->addArrow(*v_iter, newPoint, color++);
    }
    for(sieve_type::traits::coneSequence::iterator v_iter = fBegin; v_iter != fEnd;
        ++v_iter) {
      if (debug)
        std::cout << "    vertex " << vertexRenumber[*v_iter] << std::endl;
      sieve->addArrow(vertexRenumber[*v_iter], newPoint, color++);
    }
    mesh->setValue(material, newPoint, materialId);
  } // for
  mesh->stratify();
  if (debug)
    mesh->view("Mesh with Cohesive Elements");

  // Fix coordinates
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  const ALE::Obj<Mesh::label_sequence>& fVertices2 = fault->depthStratum(0);

  for(Mesh::label_sequence::iterator v_iter = fVertices2->begin();
      v_iter != fVertices2->end();
      ++v_iter) {
    coordinates->addPoint(vertexRenumber[*v_iter],
			  coordinates->getFiberDimension(*v_iter));
  } // for
  mesh->reallocate(coordinates);
  for(Mesh::label_sequence::iterator v_iter = fVertices2->begin();
      v_iter != fVertices2->end();
      ++v_iter)
    coordinates->updatePoint(vertexRenumber[*v_iter], 
			     coordinates->restrictPoint(*v_iter));

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


// End of file
