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

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc
#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array
#include <petscmat.h> // USES PETSc Mat

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES CoordSys

#include <Distribution.hh> // USES completeSection
#include <Selection.hh> // Algorithms for submeshes

#include <math.h> // USES pow(), sqrt()
#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _eqsrc(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  _eqsrc = 0; // Don't manage memory for eq source
} // destructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrc(EqKinSrc* src)
{ // eqsrc
  _eqsrc = src; // Don't manage memory for eq source
} // eqsrc

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const ALE::Obj<ALE::Mesh>& mesh,
					     const spatialdata::geocoords::CoordSys* cs,
					     const double_array& upDir,
					     const double_array& normalDir,
					     spatialdata::spatialdb::SpatialDB* matDB)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _eqsrc);
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");
  if (3 != normalDir.size())
    throw std::runtime_error("Normal direction for fault orientation must be "
			     "a vector with 3 components.");

  CohesiveTopology::createParallel(&_faultMesh, &_cohesiveToFault, mesh, id(),
				   _useLagrangeConstraints());

  //_faultMesh->view("FAULT MESH");

  // Establish pairing between constraint vertices and first cell they
  // appear in to prevent overlap in integrating Jacobian
  _calcVertexCellPairs();

  // Setup pseudo-stiffness of cohesive cells to improve conditioning
  // of Jacobian matrix
  _calcConditioning(cs, matDB);

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir, normalDir);

  _eqsrc->initialize(_faultMesh, cs);
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
				const ALE::Obj<real_section_type>& residual,
				const double t,
				topology::FieldsManager* const fields,
				const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make 2 contributions to the residual:
  //
  //   * DOF i and j: internal forces in soln field associated with 
  //                  slip
  //   * DOF k: slip values

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim*spaceDim;
  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3*numConstraintVert; // cohesive cell
  int_array cellConstraintCell(numConstraintVert);
  double_array cellOrientation(numConstraintVert*orientationSize);
  double_array cellResidual(numCorners*spaceDim);
  double_array cellSoln(numCorners*spaceDim);
  double_array cellSlip(numConstraintVert*spaceDim);
  double_array cellStiffness(numConstraintVert);

  // Get cohesive cells
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get section information
  const ALE::Obj<real_section_type>& solution = fields->getSolution();
  assert(!solution.isNull());  

  if (!_useSolnIncr) {
    // Compute slip field at current time step
    assert(0 != _eqsrc);
    _slip = _eqsrc->slip(t, _faultMesh);
    assert(!_slip.isNull());
  } else {
    // Compute increment of slip field at current time step
    assert(0 != _eqsrc);
    _slip = _eqsrc->slipIncr(t-_dt, t, _faultMesh);
    assert(!_slip.isNull());
  } // else
  
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const Mesh::point_type c_fault = _cohesiveToFault[*c_iter];

    cellResidual = 0.0;
    // Get Lagrange constraint / fault cell pairings.
    _faultMesh->restrict(_faultVertexCell, c_fault, &cellConstraintCell[0], 
			 cellConstraintCell.size());
    
    // Get orientations at fault cell's vertices.
    _faultMesh->restrict(_orientation, c_fault, &cellOrientation[0], 
		   cellOrientation.size());
    
    // Get pseudo stiffness at fault cell's vertices.
    _faultMesh->restrict(_pseudoStiffness, c_fault, &cellStiffness[0], 
			 cellStiffness.size());
    
    // Get slip at fault cell's vertices.
    _faultMesh->restrict(_slip, c_fault, &cellSlip[0], cellSlip.size());

    // Get solution at cohesive cell's vertices.
    mesh->restrict(solution, *c_iter, &cellSoln[0], cellSoln.size());
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Skip setting values if they are set by another cell
      if (cellConstraintCell[iConstraint] != c_fault)
	continue;
      
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      const double pseudoStiffness = cellStiffness[iConstraint];

      // Set slip values in residual vector; contributions are at DOF of
      // constraint vertices (k) of the cohesive cells
      for (int iDim=0; iDim < spaceDim; ++iDim)
	cellResidual[indexK*spaceDim+iDim] = 
	  cellSlip[iConstraint*spaceDim+iDim];
      
      if (_useSolnIncr) {
	// Get orientation at constraint vertex
	const real_section_type::value_type* constraintOrient = 
	  &cellOrientation[iConstraint*orientationSize];
	assert(0 != constraintOrient);

	// Entries associated with constraint forces applied at node i
	for (int iDim=0; iDim < spaceDim; ++iDim)
	  for (int kDim=0; kDim < spaceDim; ++kDim)
	    cellResidual[indexI*spaceDim+iDim] -=
	      cellSoln[indexK*spaceDim+kDim] * 
	      -constraintOrient[kDim*spaceDim+iDim] * pseudoStiffness;
	
	// Entries associated with constraint forces applied at node j
	for (int jDim=0; jDim < spaceDim; ++jDim)
	  for (int kDim=0; kDim < spaceDim; ++kDim)
	    cellResidual[indexJ*spaceDim+jDim] -=
	      cellSoln[indexK*spaceDim+kDim] * 
	      constraintOrient[kDim*spaceDim+jDim] * pseudoStiffness;
      } // if
    } // for

#if 0
    std::cout << "Updating fault residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < numConstraintVert; ++i) {
      std::cout << "  stif["<<i<<"]: " << cellStiffness[i] << std::endl;
    }
    for(int i = 0; i < numConstraintVert*spaceDim; ++i) {
      std::cout << "  slip["<<i<<"]: " << cellSlip[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  soln["<<i<<"]: " << cellSoln[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  v["<<i<<"]: " << cellResidual[i] << std::endl;
    }
#endif

    // Assemble cell contribution into field
    mesh->updateAdd(residual, *c_iter, &cellResidual[0]);
  } // for
  PetscLogFlopsNoCheck(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*7);
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveKin::integrateJacobian(
				    PetscMat* mat,
				    const double t,
				    topology::FieldsManager* const fields,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != mat);
  assert(0 != fields);
  assert(!mesh.isNull());

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  PetscErrorCode err = 0;

  // Get cohesive cells
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get section information
  const ALE::Obj<real_section_type>& solution = fields->getSolution();
  assert(!solution.isNull());  

  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim*spaceDim;

  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3*numConstraintVert; // cohesive cell
  double_array cellMatrix(numCorners*spaceDim * numCorners*spaceDim);
  double_array cellOrientation(numConstraintVert*orientationSize);
  int_array cellConstraintCell(numConstraintVert);
  double_array cellStiffness(numConstraintVert);

  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const Mesh::point_type c_fault = _cohesiveToFault[*c_iter];

    cellMatrix = 0.0;
    // Get Lagrange constraint / fault cell pairings.
    _faultMesh->restrict(_faultVertexCell, c_fault, &cellConstraintCell[0], 
			 cellConstraintCell.size());

    // Get orientations at fault cell's vertices.
    _faultMesh->restrict(_orientation, c_fault, &cellOrientation[0], 
			 cellOrientation.size());

    // Get pseudo stiffness at fault cell's vertices.
    _faultMesh->restrict(_pseudoStiffness, c_fault, &cellStiffness[0], 
			 cellStiffness.size());
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Skip setting values if they are set by another cell
      if (cellConstraintCell[iConstraint] != c_fault)
	continue;
      
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      // Get orientation at constraint vertex
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[iConstraint*orientationSize];
      assert(0 != constraintOrient);

      const double pseudoStiffness = cellStiffness[iConstraint];

      // Scale orientation information by pseudo-stiffness to bring
      // constraint forces in solution vector to the same order of
      // magnitude as the displacements to prevent ill-conditioning

      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = indexI*spaceDim+iDim;
	  const int col = indexK*spaceDim+kDim;
	  cellMatrix[row*numCorners*spaceDim+col] =
	    -constraintOrient[kDim*spaceDim+iDim]*pseudoStiffness;
	  cellMatrix[col*numCorners*spaceDim+row] =
	    -constraintOrient[kDim*spaceDim+iDim];
	} // for

      // Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = indexJ*spaceDim+jDim;
	  const int col = indexK*spaceDim+kDim;
	  cellMatrix[row*numCorners*spaceDim+col] =
	    constraintOrient[kDim*spaceDim+jDim]*pseudoStiffness;
	  cellMatrix[col*numCorners*spaceDim+row] =
	    constraintOrient[kDim*spaceDim+jDim];
	} // for
    } // for

    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", solution);
    // Note: We are not really adding values because we prevent
    // overlap across cells. We use ADD_VALUES for compatibility with
    // the other integrators.
    err = updateOperator(*mat, mesh, solution, globalOrder,
			 *c_iter, &cellMatrix[0], ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for
  PetscLogFlopsNoCheck(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*4);
  _needNewJacobian = false;
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveKin::verifyConfiguration(
					      const ALE::Obj<Mesh>& mesh)
{ // verifyConfiguration
  assert(0 != _quadrature);

  // check compatibility of mesh and quadrature scheme
  const int dimension = mesh->getDimension()-1;
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Dimension of reference cell in quadrature scheme (" 
	<< _quadrature->cellDim() 
	<< ") does not match dimension of cells in mesh (" 
	<< dimension << ") for fault '" << label()
	<< "'.";
    throw std::runtime_error(msg.str());
  } // if

  const int numCorners = _quadrature->refGeometry().numCorners();
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    mesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  for (Mesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieve->nCone(*c_iter, mesh->depth())->size();
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
const ALE::Obj<pylith::real_section_type>&
pylith::faults::FaultCohesiveKin::vertexField(
				    VectorFieldEnum* fieldType,
				    const char* name,
				    const ALE::Obj<Mesh>& mesh,
				    topology::FieldsManager* fields)
{ // vertexField
  assert(!_faultMesh.isNull());
  assert(!_orientation.isNull());
  assert(0 != _eqsrc);

  const int cohesiveDim = _faultMesh->getDimension();

  if (0 == strcasecmp("slip", name)) {
    *fieldType = VECTOR_FIELD;
    return _slip;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    _bufferVertexVector = _orientation->getFibration(0);
    *fieldType = VECTOR_FIELD;
    return _bufferVertexVector;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    _bufferVertexVector = _orientation->getFibration(1);
    *fieldType = VECTOR_FIELD;
    return _bufferVertexVector;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const int space = 
      (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    _bufferVertexVector = _orientation->getFibration(space);
    *fieldType = VECTOR_FIELD;
    return _bufferVertexVector;

  } else if (0 == strcasecmp("final_slip", name)) {
    _bufferVertexVector = _eqsrc->finalSlip();
    *fieldType = VECTOR_FIELD;
    return _bufferVertexVector;
  } else if (0 == strcasecmp("slip_time", name)) {
    _bufferVertexScalar = _eqsrc->slipTime();
    *fieldType = SCALAR_FIELD;
    return _bufferVertexScalar;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name
	<< "' for fault '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Return generic section to satisfy member function definition.
  return _bufferVertexScalar;
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::FaultCohesiveKin::cellField(
				    VectorFieldEnum* fieldType,
				    const char* name,
				    const ALE::Obj<Mesh>& mesh,
				    topology::FieldsManager* fields)
{ // cellField
  assert(!_faultMesh.isNull());
  assert(!_orientation.isNull());
  assert(0 != fields);
  assert(0 != _eqsrc);

  const int cohesiveDim = _faultMesh->getDimension();

  if (0 == strcasecmp("traction_change", name)) {
    _allocateBufferCellVector();
    *fieldType = VECTOR_FIELD;
    const ALE::Obj<real_section_type>& solution = fields->getSolution();
    _calcTractionsChange(&_bufferCellVector, solution);
    return _bufferCellVector;
  } // if

  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown vertex field '" << name
      << "' for fault '" << label() << "'.";
  throw std::runtime_error(msg.str());

  // Return generic section to satisfy member function definition.
  return _bufferCellVector;
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveKin::_calcOrientation(const double_array& upDir,
						   const double_array& normalDir)
{ // _calcOrientation
  assert(!_faultMesh.isNull());

  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  // Create orientation section for fault (constraint) vertices
  const int cohesiveDim = _faultMesh->getDimension();
  const int spaceDim = cohesiveDim + 1;
  const int orientationSize = spaceDim*spaceDim;
  _orientation = new real_section_type(_faultMesh->comm(), 
				       _faultMesh->debug());
  assert(!_orientation.isNull());
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    _orientation->addSpace();
  assert(cohesiveDim+1 == _orientation->getNumSpaces());
  _orientation->setFiberDimension(vertices, orientationSize);
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    _orientation->setFiberDimension(vertices, spaceDim, iDim);
  _faultMesh->allocate(_orientation);
  
  // Compute orientation of fault at constraint vertices

  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    _faultMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Set orientation function
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());

  // Loop over cohesive cells, computing orientation at constraint vertices
  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array vertexOrientation(orientationSize);
  double_array faceVertices(numBasis*spaceDim);
  
  // Get fault cells (1 dimension lower than top-level cells)
  const ALE::Obj<Mesh::label_sequence>& cells = 
    _faultMesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<sieve_type>& sieve = _faultMesh->getSieve();
  assert(!sieve.isNull());
  const int faultDepth = _faultMesh->depth();  // depth of fault cells
  typedef ALE::SieveAlg<Mesh> SieveAlg;

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    _faultMesh->restrict(coordinates, *c_iter, 
			 &faceVertices[0], faceVertices.size());

    const ALE::Obj<SieveAlg::coneArray>& cone =
      SieveAlg::nCone(_faultMesh, *c_iter, faultDepth);
    assert(!cone.isNull());
    const SieveAlg::coneArray::iterator vBegin = cone->begin();
    const SieveAlg::coneArray::iterator vEnd = cone->end();

    int iBasis = 0;
    for(SieveAlg::coneArray::iterator v_iter=vBegin;
	v_iter != vEnd;
	++v_iter, ++iBasis) {
      // Compute Jacobian and determinant of Jacobian at vertex
      double_array vertex(&verticesRef[iBasis*cohesiveDim], cohesiveDim);
      cellGeometry.jacobian(&jacobian, &jacobianDet, faceVertices, vertex);

      // Compute orientation
      cellGeometry.orientation(&vertexOrientation, jacobian, jacobianDet, 
			       upDir);
      
      // Update orientation
      _orientation->updateAddPoint(*v_iter, &vertexOrientation[0]);
    } // for
  } // for

  // Assemble orientation information
  ALE::Distribution<Mesh>::completeSection(_faultMesh, _orientation);

  // Loop over vertices, make orientation information unit magnitude
  double_array vertexDir(orientationSize);
  int count = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    const real_section_type::value_type* vertexOrient = 
      _orientation->restrictPoint(*v_iter);
    assert(0 != vertexOrient);

    assert(spaceDim*spaceDim == orientationSize);
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      double mag = 0;
      for (int jDim=0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	mag += pow(vertexOrient[index+jDim],2);
      mag = sqrt(mag);
      assert(mag > 0.0);
      for (int jDim=0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	vertexDir[index+jDim] = 
	  vertexOrient[index+jDim] / mag;
    } // for

    _orientation->updatePoint(*v_iter, &vertexDir[0]);
  } // for
  PetscLogFlopsNoCheck(count * orientationSize * 4);

  if (2 == cohesiveDim) {
    // Check orientation of first vertex, if dot product of fault
    // normal with preferred normal is negative, flip up/down dip direction.
    // If the user gives the correct normal direction, we should end
    // up with left-lateral-slip, reverse-slip, and fault-opening for
    // positive slip values.

    const real_section_type::value_type* vertexOrient = 
      _orientation->restrictPoint(*vertices->begin());
    assert(0 != vertexOrient);

    double_array vertNormalDir(&vertexOrient[6], 3);
    const double dot = 
      normalDir[0]*vertNormalDir[0] +
      normalDir[1]*vertNormalDir[1] +
      normalDir[2]*vertNormalDir[2];
    if (dot < 0.0)
      for (Mesh::label_sequence::iterator v_iter=vertices->begin();
	   v_iter != verticesEnd;
	   ++v_iter) {
	const real_section_type::value_type* vertexOrient = 
	  _orientation->restrictPoint(*v_iter);
	assert(0 != vertexOrient);
	assert(9 == _orientation->getFiberDimension(*v_iter));
	// Keep along-strike direction
	for (int iDim=0; iDim < 3; ++iDim)
	  vertexDir[iDim] = vertexOrient[iDim];
	// Flip up-dip direction
	for (int iDim=3; iDim < 6; ++iDim)
	  vertexDir[iDim] = -vertexOrient[iDim];
	// Keep normal direction
	for (int iDim=6; iDim < 9; ++iDim)
	  vertexDir[iDim] = vertexOrient[iDim];
	
	// Update direction
	_orientation->updatePoint(*v_iter, &vertexDir[0]);
      } // for

    PetscLogFlopsNoCheck(5 + count * 3);
  } // if

  //_orientation->view("ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
// Calculate conditioning field.
void
pylith::faults::FaultCohesiveKin::_calcVertexCellPairs(void)
{ // _calcVertexCellPairs
  assert(!_faultMesh.isNull());

  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  
  const int noCell = -1;
  _faultVertexCell = new int_section_type(_faultMesh->comm(), 
					  _faultMesh->debug());
  assert(!_faultVertexCell.isNull());
  _faultVertexCell->setFiberDimension(vertices, 1);
  _faultMesh->allocate(_faultVertexCell);
  
  // Set values to noCell
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter)
    _faultVertexCell->updatePoint(*v_iter, &noCell);
  
  const ALE::Obj<Mesh::label_sequence>& cells = 
    _faultMesh->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<sieve_type>& sieve = _faultMesh->getSieve();
  assert(!sieve.isNull());
  const int faultDepth = _faultMesh->depth();  // depth of fault cells
  typedef ALE::SieveAlg<Mesh> SieveAlg;
  
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const ALE::Obj<SieveAlg::coneArray>& cone =
      SieveAlg::nCone(_faultMesh, *c_iter, faultDepth);
    assert(!cone.isNull());
    const SieveAlg::coneArray::iterator vBegin = cone->begin();
    const SieveAlg::coneArray::iterator vEnd = cone->end();
    const int coneSize = cone->size();

    // If haven't set cell-constraint pair, then set it for current
    // cell, otherwise move on.
    SieveAlg::coneArray::iterator v_iter = vBegin;
    for(int i=0; i < coneSize; ++i, ++v_iter) {
      const int_section_type::value_type* curCell = 
	_faultVertexCell->restrictPoint(*v_iter);
      assert(0 != curCell);
      if (noCell == *curCell) {
	int point = *c_iter;
	_faultVertexCell->updatePoint(*v_iter, &point);
      } // if
    } // for
  } // for

  //_faultVertexCell->view("VERTEX/CELL PAIRINGS");
} // _calcVertexCallPairs

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveKin::_calcConditioning(
				 const spatialdata::geocoords::CoordSys* cs,
				 spatialdata::spatialdb::SpatialDB* matDB)
{ // _calcConditioning
  assert(0 != cs);
  assert(0 != matDB);
  assert(!_faultMesh.isNull());

  const int spaceDim = cs->spaceDim();

  // Get vertices in fault mesh
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  
  _pseudoStiffness = new real_section_type(_faultMesh->comm(), 
					   _faultMesh->debug());
  assert(!_pseudoStiffness.isNull());
  _pseudoStiffness->setFiberDimension(vertices, 1);
  _faultMesh->allocate(_pseudoStiffness);
  
  matDB->open();
  const char* stiffnessVals[] = { "density", "vs" };
  const int numStiffnessVals = 2;
  matDB->queryVals(stiffnessVals, numStiffnessVals);
  
  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    _faultMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  double_array matprops(numStiffnessVals);
  int count = 0;
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    const real_section_type::value_type* vertexCoords = 
      coordinates->restrictPoint(*v_iter);
    assert(0 != vertexCoords);
    int err = matDB->query(&matprops[0], numStiffnessVals, vertexCoords, 
			   spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find material properties at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vertexCoords[i];
      msg << ") using spatial database " << matDB->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    
    const double density = matprops[0];
    const double vs = matprops[1];
    const double mu = density * vs*vs;
    //const double mu = 1.0;
    _pseudoStiffness->updatePoint(*v_iter, &mu);
  } // for
  PetscLogFlopsNoCheck(count * 2);
} // _calcConditioning

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveKin::_calcTractionsChange(
				 ALE::Obj<real_section_type>* tractions,
				 const ALE::Obj<real_section_type>& solution)
{ // _calcTractionsChange
  assert(0 != tractions);
  assert(!tractions->isNull());
  assert(!solution.isNull());
  assert(!_faultMesh.isNull());
  assert(!_pseudoStiffness.isNull());
  assert(0 != _quadrature);

  const ALE::Obj<Mesh::label_sequence>& cells = 
    _faultMesh->heightStratum(0);
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& basis = _quadrature->basis();

  double_array solutionCell(numBasis*spaceDim);
  double_array stiffnessCell(numBasis);
  double_array tractionsCell(numQuadPts*spaceDim);

  int count = 0;
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter, ++count) {
    _faultMesh->restrict(solution, *c_iter, 
			 &solutionCell[0], solutionCell.size());
    _faultMesh->restrict(_pseudoStiffness, *c_iter, 
			 &stiffnessCell[0], stiffnessCell.size());
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int iBasis=0; iBasis < numBasis; ++iBasis)
	for (int iDim=0; iDim < spaceDim; ++iDim)
	  tractionsCell[iQuad*spaceDim+iDim] = basis[iBasis] *
	    solutionCell[iBasis*spaceDim+iDim] * stiffnessCell[iBasis];
    (*tractions)->updatePoint(*c_iter, &tractionsCell[0]);
  } // for

  //solution->view("SOLUTION");
  //(*tractions)->view("TRACTIONS");
} // _calcTractionsChange

// ----------------------------------------------------------------------
// Allocate scalar field for output of vertex information.
void
pylith::faults::FaultCohesiveKin::_allocateBufferVertexScalar(void)
{ // _allocateBufferVertexScalar
  const int fiberDim = 1;
  if (_bufferVertexScalar.isNull()) {
    _bufferVertexScalar = new real_section_type(_faultMesh->comm(), 
						_faultMesh->debug());
    const ALE::Obj<Mesh::label_sequence>& vertices = 
      _faultMesh->depthStratum(0);
    _bufferVertexScalar->setFiberDimension(vertices, fiberDim);
    _faultMesh->allocate(_bufferVertexScalar);
  } // if
} // _allocateBufferVertexScalar

// ----------------------------------------------------------------------
// Allocate vector field for output of vertex information.
void
pylith::faults::FaultCohesiveKin::_allocateBufferVertexVector(void)
{ // _allocateBufferVertexVector
  assert(0 != _quadrature);
  const int fiberDim = _quadrature->spaceDim();
  if (_bufferVertexVector.isNull()) {
    _bufferVertexVector = new real_section_type(_faultMesh->comm(), 
						_faultMesh->debug());
    const ALE::Obj<Mesh::label_sequence>& vertices = 
      _faultMesh->depthStratum(0);
    _bufferVertexVector->setFiberDimension(vertices, fiberDim);
    _faultMesh->allocate(_bufferVertexVector);
  } // if  
} // _allocateBufferVertexVector

// ----------------------------------------------------------------------
// Allocate vector field for output of cell information.
void
pylith::faults::FaultCohesiveKin::_allocateBufferCellVector(void)
{ // _allocateBufferCellVector
  assert(0 != _quadrature);
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const int fiberDim = numQuadPts * spaceDim;
  if (_bufferCellVector.isNull()) {
    _bufferCellVector = new real_section_type(_faultMesh->comm(), 
					      _faultMesh->debug());
    const ALE::Obj<Mesh::label_sequence>& cells = 
      _faultMesh->heightStratum(0);
    _bufferCellVector->setFiberDimension(cells, fiberDim);
    _faultMesh->allocate(_bufferCellVector);
  } // if  
} // _allocateBufferCellVector


// End of file 
