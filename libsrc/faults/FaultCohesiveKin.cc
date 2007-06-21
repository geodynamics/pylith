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

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager
#include "pylith/utils/array.hh" // USES double_array
#include <petscmat.h> // USES PETSc Mat

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES CoordSys

#include <Distribution.hh> // USES completeSection

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

  /* First find vertices associated with Lagrange multiplier
   * constraints in cohesive cells and compute the orientation of the
   * fault at these locations.
   */

  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  // Create set of vertices associated with Lagrange multiplier constraints
  _constraintVert.clear();
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    // Vertices for each cohesive cell are in three groups of N.
    // 0 to N-1: vertices on negative side of the fault
    // N-1 to 2N-1: vertices on positive side of the fault
    // 2N to 3N-1: vertices associated with constraint forces
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    assert(!cone.isNull());
    const sieve_type::traits::coneSequence::iterator vBegin = cone->begin();
    const sieve_type::traits::coneSequence::iterator vEnd = cone->end();
    const int coneSize = cone->size();
    assert(coneSize % 3 == 0);
    sieve_type::traits::coneSequence::iterator v_iter = vBegin;
    // Skip over non-constraint vertices
    for (int i=0, numSkip=2*coneSize/3; i < numSkip; ++i)
      ++v_iter;
    // Add constraint vertices to set
    for(int i=0, numConstraintVert=coneSize/3; 
	i < numConstraintVert; 
	++i, ++v_iter)
      _constraintVert.insert(*v_iter);
  } // for

  // Create orientation section for constraint vertices
  const int spaceDim = cs->spaceDim();
  const int orientationSize = spaceDim*spaceDim;
  _orientation = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_orientation.isNull());
  const std::set<Mesh::point_type>::const_iterator vertConstraintBegin = 
    _constraintVert.begin();
  const std::set<Mesh::point_type>::const_iterator vertConstraintEnd = 
    _constraintVert.end();
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter)
    _orientation->setFiberDimension(*v_iter, orientationSize);
  mesh->allocate(_orientation);

  // Compute orientation of fault at constraint vertices

  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Set orientation function
  const int cohesiveDim = mesh->getDimension()-1;
  assert(cohesiveDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());
  orient_fn_type orientFn;
  switch (cohesiveDim)
    { // switch
    case 0 :
      orientFn = _orient1D;
      break;
    case 1 :
      orientFn = _orient2D;
      break;
    case 2 :
      orientFn = _orient3D;
      break;
    default :
      assert(0);
    } // switch

  // Loop over cohesive cells, computing orientation at constraint vertices
  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = _quadrature->vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array vertexOrientation(orientationSize);
  double_array cohesiveVertices(3*numBasis*spaceDim);
  double_array faceVertices(numBasis*spaceDim);

  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    mesh->restrict(coordinates, *c_iter, 
		   &cohesiveVertices[0], cohesiveVertices.size());
    const int size = numBasis*spaceDim;
    for (int i=0, offset=2*size; i < size; ++i)
      faceVertices[i] = cohesiveVertices[offset+i];

    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    assert(!cone.isNull());
    const sieve_type::traits::coneSequence::iterator vBegin = cone->begin();
    const sieve_type::traits::coneSequence::iterator vEnd = cone->end();

    // skip over non-constraint vertices
    sieve_type::traits::coneSequence::iterator vConstraintBegin = vBegin;
    const int numSkip = 2*numBasis;
    for (int i=0; i < numSkip; ++i)
      ++vConstraintBegin;

    int iBasis = 0;
    for(sieve_type::traits::coneSequence::iterator v_iter=vConstraintBegin;
	v_iter != vEnd;
	++v_iter, ++iBasis) {
      // Compute Jacobian and determinant of Jacobian at vertex
      double_array vertex(&verticesRef[iBasis*cohesiveDim], cohesiveDim);
      cellGeometry.jacobian(&jacobian, &jacobianDet, faceVertices, vertex);

      // Compute orientation
      orientFn(&vertexOrientation, jacobian, jacobianDet, upDir);
      
      // Update orientation
      _orientation->updateAddPoint(*v_iter, &vertexOrientation[0]);
    } // for
  } // for

  // Assemble orientation information
  ALE::Distribution<Mesh>::completeSection(mesh, _orientation);

  // Loop over vertices, make orientation information unit magnitude
  double_array vertexDir(orientationSize);
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter) {
    const real_section_type::value_type* vertexOrient = 
      _orientation->restrictPoint(*v_iter);
    
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

  if (2 == cohesiveDim) {
    // Check orientation of first vertex, if dot product of fault
    // normal with preferred normal is negative, flip all along-strike
    // and fault-normal directions (keep up-dip direction the same).
    // If the user gives the correct normal direction, we should end
    // up with left-lateral-slip, reverse-slip, and fault-opening for
    // positive slip values.

    const real_section_type::value_type* vertexOrient = 
      _orientation->restrictPoint(*vertConstraintBegin);
    assert(9 == _orientation->getFiberDimension(*vertConstraintBegin));
    double_array vertNormalDir(&vertexOrient[6], 3);
    const double dot = 
      normalDir[0]*vertNormalDir[0] +
      normalDir[1]*vertNormalDir[1] +
      normalDir[2]*vertNormalDir[2];
    if (dot < 0.0)
      for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
	   v_iter != vertConstraintEnd;
	   ++v_iter) {
	const real_section_type::value_type* vertexOrient = 
	  _orientation->restrictPoint(*v_iter);
	assert(9 == _orientation->getFiberDimension(*v_iter));
	// Flip along-strike direction
	for (int iDim=0; iDim < 3; ++iDim)
	  vertexDir[iDim] = -vertexOrient[iDim];
	// Keep up-dip direction
	for (int iDim=3; iDim < 6; ++iDim)
	  vertexDir[iDim] = vertexOrient[iDim];
	// Flip normal direction
	for (int iDim=6; iDim < 9; ++iDim)
	  vertexDir[iDim] = -vertexOrient[iDim];
	
	// Update direction
	_orientation->updatePoint(*v_iter, &vertexDir[0]);
      } // for
  } // if

  _eqsrc->initialize(mesh, *_faultMesh, _constraintVert, cs);

  // Establish pairing between constraint vertices and first cell they
  // appear in to prevent overlap in integrating Jacobian
  const int noCell = -1;
  _constraintCell = new int_section_type(mesh->comm(), mesh->debug());
  assert(!_constraintCell.isNull());
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter)
    _constraintCell->setFiberDimension(*v_iter, 1);
  mesh->allocate(_constraintCell);
  // Set values to noCell
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter)
    _constraintCell->updatePoint(*v_iter, &noCell);

  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    assert(!cone.isNull());
    const sieve_type::traits::coneSequence::iterator vBegin = cone->begin();
    const sieve_type::traits::coneSequence::iterator vEnd = cone->end();
    const int coneSize = cone->size();
    assert(coneSize % 3 == 0);
    sieve_type::traits::coneSequence::iterator v_iter = vBegin;
    // Skip over non-constraint vertices
    for (int i=0, numSkip=2*coneSize/3; i < numSkip; ++i)
      ++v_iter;
    // If haven't set cell-constraint pair, then set it for current
    // cell, otherwise move on.
    for(int i=0, numConstraintVert=coneSize/3; 
	i < numConstraintVert; 
	++i, ++v_iter) {
      const int_section_type::value_type* curCell = 
	_constraintCell->restrictPoint(*v_iter);
      if (noCell == *curCell) {
	int point = *c_iter;
	_constraintCell->updatePoint(*v_iter, &point);
      } // if
    } // for
  } // for

  // Setup pseudo-stiffness of cohesive cells to improve conditioning
  // of Jacobian matrix
  _pseudoStiffness = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_pseudoStiffness.isNull());
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter)
    _pseudoStiffness->setFiberDimension(*v_iter, 1);
  mesh->allocate(_pseudoStiffness);
    
  matDB->open();
  const char* stiffnessVals[] = { "density", "vs" };
  const int numStiffnessVals = 2;
  matDB->queryVals(stiffnessVals, numStiffnessVals);
  
  double_array matprops(numStiffnessVals);
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertConstraintBegin;
       v_iter != vertConstraintEnd;
       ++v_iter) {
    const real_section_type::value_type* vertexCoords = 
      coordinates->restrictPoint(*v_iter);
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
    _pseudoStiffness->updatePoint(*v_iter, &mu);
  } // for
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

  // Get section information
  const ALE::Obj<real_section_type>& solution = fields->getSolution();
  assert(!solution.isNull());  

  ALE::Obj<real_section_type> slip;
  if (!_useSolnIncr) {
    // Compute slip field at current time step
    assert(0 != _eqsrc);
    slip = _eqsrc->slip(t, _constraintVert);
    assert(!slip.isNull());
  } else {
    // Compute increment of slip field at current time step
    assert(0 != _eqsrc);
    slip = _eqsrc->slipIncr(t-_dt, t, _constraintVert);
    assert(!slip.isNull());
  } // else
  
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    cellResidual = 0.0;
    // Get constraint/cell pairings (only valid at constraint vertices)
    mesh->restrict(_constraintCell, *c_iter, &cellConstraintCell[0], 
		   cellConstraintCell.size());
    
    // Get orientations at cells vertices (only valid at constraint vertices)
    mesh->restrict(_orientation, *c_iter, &cellOrientation[0], 
		   cellOrientation.size());
    
    // Get pseudo stiffness at cells vertices (only valid at
    // constraint vertices)
    mesh->restrict(_pseudoStiffness, *c_iter, &cellStiffness[0], 
		   cellStiffness.size());
    
    // Get slip at cells vertices (only valid at constraint vertices)
    mesh->restrict(slip, *c_iter, &cellSlip[0], cellSlip.size());

    // Get slip at cells vertices (valid at all cohesive vertices)
    mesh->restrict(solution, *c_iter, &cellSoln[0], cellSoln.size());
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Skip setting values if they are set by another cell
      if (cellConstraintCell[iConstraint] != *c_iter)
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
	  cellSlip[iConstraint*spaceDim+iDim] * pseudoStiffness;
      
      // Get orientation at constraint vertex
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[iConstraint*orientationSize];

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
    } // for

    // Assemble cell contribution into field
    mesh->updateAdd(residual, *c_iter, &cellResidual[0]);
  } // for
  PetscErrorCode err = PetscLogFlops(numConstraintVert*spaceDim*spaceDim*6 +
				     numConstraintVert*spaceDim*2);
  if (err)
    throw std::runtime_error("Logging PETSc flops failed.");
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
    cellMatrix = 0.0;
    // Get orientations at cells vertices (only valid at constraint vertices)
    mesh->restrict(_orientation, *c_iter, &cellOrientation[0], 
		   cellOrientation.size());

    // Get pseudo stiffness at cells vertices (only valid at
    // constraint vertices)
    mesh->restrict(_pseudoStiffness, *c_iter, &cellStiffness[0], 
		   cellStiffness.size());
    
    // Get constraint/cell pairings (only valid at constraint vertices)
    mesh->restrict(_constraintCell, *c_iter, &cellConstraintCell[0], 
		   cellConstraintCell.size());

    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Skip setting values if they are set by another cell
      if (cellConstraintCell[iConstraint] != *c_iter)
	continue;
      
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      // Get orientation at constraint vertex
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[iConstraint*orientationSize];

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
	    cellMatrix[row*numCorners*spaceDim+col]; // symmetric
	} // for

      // Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = indexJ*spaceDim+jDim;
	  const int col = indexK*spaceDim+kDim;
	  cellMatrix[row*numCorners*spaceDim+col] =
	    constraintOrient[kDim*spaceDim+jDim]*pseudoStiffness;
	  cellMatrix[col*numCorners*spaceDim+row] =
	    cellMatrix[row*numCorners*spaceDim+col]; // symmetric
	} // for
    } // for
    err = PetscLogFlops(numConstraintVert*spaceDim*spaceDim*4);
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

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
  _needNewJacobian = false;
} // integrateJacobian
  

// End of file 
