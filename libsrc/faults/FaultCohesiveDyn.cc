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

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::FaultCohesiveDyn::deallocate(void)
{ // deallocate
  FaultCohesive::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					     const double upDir[3],
					     const double normalDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != normalDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  // _calcOrientation()
  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir, normalDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

  //_getInitialTractions()
  
  
  // _initializeConstitutiveModel()
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(
			   const topology::Field<topology::Mesh>& residual,
			   const double t,
			   topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _faultMesh); // changed to fault mesh
  assert(0 != _fields); // change to fields

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numCornersCohesive = 2*numBasis;

  // Allocate vectors for cell values.
  double_array tractionsCell(numQuadPts*spaceDim);
  double_array residualCell(numCornersCohesive*spaceDim);

  // GET COHESIVE CELLS (see FaultCohesiveKin)
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  // Add setting cellsBegin

  // Get cell information
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh(); // updated (subSieveMesh -> faultSieveMesh)
  const SieveSubMesh::label_sequence::iterator cellsEnd = cellsCohesive->end();
  assert(!faultSieveMesh.isNull());

  // GET PARAMETERS FOR FAULT CONSTITUTIVE MODEL (later)
  // TEMPORARY hardwire to (1,0,0)
  for (int i=0; i < numQuadPts; ++i)
    tractionsCell[i*spaceDim] = 1.0;

  // Get sections (see FaultCohesiveKin)
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						      &residualCell[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // CHANGE TO LOOP OVER COHESIVE CELLS (See FaultCohesiveKin)
 
  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsCohesive->begin();
       c_iter != cellsEnd; // cellsEnd
       ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    residualCell = 0.0;

#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

    // USE FAULT CONSTITUTIVE MODEL TO GET TRACTION (later)
    
    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            residualCell[iBasis*spaceDim+iDim] += 
	      tractionsCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateAdd(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveDyn::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const double t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  // See Neumann
  throw std::logic_error("FaultCohesiveDyn::integrateJacobian() not implemented.");
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveDyn::verifyConfiguration(
					    const topology::Mesh& mesh) const
{ // verifyConfiguration
  assert(0 != _quadrature);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(label())) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label()
	<< " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if  

  // check compatibility of mesh and quadrature scheme
  const int dimension = mesh.dimension()-1;
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Dimension of reference cell in quadrature scheme (" 
	<< _quadrature->cellDim() 
	<< ") does not match dimension of cells in mesh (" 
	<< dimension << ") for fault '" << label()
	<< "'.";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieveMesh->getNumCellCorners(*c_iter);
    if (3*numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Number of vertices in reference cell (" << numCorners 
	  << ") is not compatible with number of vertices (" << cellNumCorners
	  << ") in cohesive cell " << *c_iter << " for fault '"
	  << label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::vertexField(
				       const char* name,
				       const topology::SolutionFields* fields)
{ // vertexField
  throw std::logic_error("FaultCohesiveDyn::vertexField() not implemented.");
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::cellField(
				      const char* name,
				      const topology::SolutionFields* fields)
{ // cellField
  throw std::logic_error("FaultCohesiveDyn::cellField() not implemented.");
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesive::_calcOrientation(const double upDir[3],
						   const double normalDir[3])
{ // _calcOrientation
  assert(0 != upDir);
  assert(0 != normalDir);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  double_array upDirArray(upDir, 3);

  // Get vertices in fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  // Containers for orientation information.
  const int cohesiveDim = _faultMesh->dimension();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim*spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  const double_array& quadWts = _quadrature->quadWts();
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientationVertex(orientationSize);
  double_array coordinatesCell(numBasis*spaceDim);
  double_array refCoordsVertex(cohesiveDim);

  // Allocate orientation field.
  _fields->add("orientation", "orientation");
  topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");
  const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  orientation.newSection(slip, orientationSize);
  const ALE::Obj<RealSection>& orientationSection = orientation.section();
  assert(!orientationSection.isNull());
  // Create subspaces for along-strike, up-dip, and normal directions
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    orientationSection->addSpace();
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    orientationSection->setFiberDimension(vertices, spaceDim, iDim);
  orientation.allocate();
  orientation.zero();
  
  // Get fault cells.
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Compute orientation of fault at constraint vertices

  // Get section containing coordinates of vertices
  const ALE::Obj<RealSection>& coordinatesSection = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinatesSection.isNull());
  topology::Mesh::RestrictVisitor coordinatesVisitor(*coordinatesSection,
						     coordinatesCell.size(),
						     &coordinatesCell[0]);

  // Set orientation function
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices
  
  const ALE::Obj<SieveSubMesh::sieve_type>& sieve = faultSieveMesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<SieveSubMesh> SieveAlg;

  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve, (size_t) pow(sieve->getMaxConeSize(), std::max(0, faultSieveMesh->depth())));

  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Get orientations at fault cell's vertices.
    coordinatesVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordinatesVisitor);

    ncV.clear();
    ALE::ISieveTraversal<SieveSubMesh::sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int               coneSize = ncV.getSize();
    const Mesh::point_type *cone     = ncV.getPoints();
    
    for (int v=0; v < coneSize; ++v) {
      // Compute Jacobian and determinant of Jacobian at vertex
      memcpy(&refCoordsVertex[0], &verticesRef[v*cohesiveDim],
	     cohesiveDim*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell,
			    refCoordsVertex);

      // Compute orientation
      cellGeometry.orientation(&orientationVertex, jacobian, jacobianDet, 
			       upDirArray);
      
      // Update orientation
      orientationSection->updateAddPoint(cone[v], &orientationVertex[0]);
    } // for
  } // for

  //orientation.view("ORIENTATION BEFORE COMPLETE");

  // Assemble orientation information
  orientation.complete();

  // Loop over vertices, make orientation information unit magnitude
  double_array vertexDir(orientationSize);
  int count = 0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    orientationVertex = 0.0;
    orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
				      orientationVertex.size());
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      double mag = 0;
      for (int jDim=0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	mag += pow(orientationVertex[index+jDim],2);
      mag = sqrt(mag);
      assert(mag > 0.0);
      for (int jDim=0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	orientationVertex[index+jDim] /= mag;
    } // for

    orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
  } // for
  PetscLogFlops(count * orientationSize * 4);

  if (2 == cohesiveDim && vertices->size() > 0) {
    // Check orientation of first vertex, if dot product of fault
    // normal with preferred normal is negative, flip up/down dip
    // direction.
    //
    // If the user gives the correct normal direction (points from
    // footwall to ahanging wall), we should end up with
    // left-lateral-slip, reverse-slip, and fault-opening for positive
    // slip values.
    //
    // When we flip the up/down dip direction, we create a left-handed
    // strike/dip/normal coordinate system, but it gives the correct
    // sense of slip. In reality the strike/dip/normal directions that
    // are used are the opposite of what we would want, but we cannot
    // flip the fault normal direction because it is tied to how the
    // cohesive cells are created.
    
    assert(vertices->size() > 0);
    orientationSection->restrictPoint(*vertices->begin(), &orientationVertex[0],
				      orientationVertex.size());
				      
    assert(3 == spaceDim);
    double_array normalDirVertex(&orientationVertex[6], 3);
    const double normalDot = 
      normalDir[0]*normalDirVertex[0] +
      normalDir[1]*normalDirVertex[1] +
      normalDir[2]*normalDirVertex[2];
    
    const int istrike = 0;
    const int idip = 3;
    const int inormal = 6;
    if (normalDot < 0.0) {
      // Flip dip direction
      for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
	   v_iter != verticesEnd;
	   ++v_iter) {
	orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
					  orientationVertex.size());
	assert(9 == orientationSection->getFiberDimension(*v_iter));
	for (int iDim=0; iDim < 3; ++iDim) // flip dip
	  orientationVertex[idip+iDim] *= -1.0;
	
	// Update direction
	orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
      } // for
      PetscLogFlops(5 + count * 3);
    } // if
  } // if

  //orientation.view("ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesive::_calcArea(void)
{ // _calcArea
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // Containers for area information
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  double jacobianDet = 0;
  double_array areaCell(numBasis);

  // Get vertices in fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  // Allocate area field.
  _fields->add("area", "area");

  topology::Field<topology::SubMesh>& area = _fields->get("area");
  const topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  area.newSection(slip, 1);
  area.allocate();
  area.zero();
  const ALE::Obj<RealSection>& areaSection = area.section();
  assert(!areaSection.isNull());
  topology::Mesh::UpdateAddVisitor areaVisitor(*areaSection, &areaCell[0]);  
  
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Loop over cells in fault mesh, compute area
  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    areaCell = 0.0;
    
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute area
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt*basis[iQuad*numBasis+iBasis];
	areaCell[iBasis] += dArea;
      } // for
    } // for
    areaVisitor.clear();
    faultSieveMesh->updateAdd(*c_iter, areaVisitor);

    PetscLogFlops( numQuadPts*(1+numBasis*2) );
  } // for

  // Assemble area information
  area.complete();

#if 0 // DEBUGGING
  area.view("AREA");
  //_faultMesh->getSendOverlap()->view("Send fault overlap");
  //_faultMesh->getRecvOverlap()->view("Receive fault overlap");
#endif
} // _calcArea


// End of file 
