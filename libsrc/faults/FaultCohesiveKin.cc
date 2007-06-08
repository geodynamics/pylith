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

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
					     const double_array& upDir)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _eqsrc);
  assert(0 != _faultMesh);
  assert(!_faultMesh->isNull());
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");

  /* First use fault mesh to get orientation of vertices in fault mesh
   * We will transfer the orientation from the fault mesh vertices to
   * the Lagrange constraint vertices.
   */

  // Allocate section for orientation of vertices in fault mesh
  ALE::Obj<real_section_type> orientation = 
    new real_section_type((*_faultMesh)->comm(), (*_faultMesh)->debug());
  assert(!orientation.isNull());
  const int cellDim = (*_faultMesh)->getDimension();
  const int spaceDim = cs->spaceDim();
  const int orientationSize = (cellDim > 0) ? cellDim*spaceDim : 1;
  orientation->setFiberDimension((*_faultMesh)->depthStratum(0), 
				 orientationSize);
  (*_faultMesh)->allocate(orientation);
  
  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Set orientation function
  assert(cellDim == _quadrature->cellDim());
  assert(spaceDim == _quadrature->spaceDim());
  orient_fn_type orientFn;
  switch (cellDim)
    { // switch
    case 1 :
      orientFn = _orient1D;
      break;
    case 2 :
      orientFn = _orient2D;
      break;
    case 3 :
      orientFn = _orient3D;
      break;
    default :
      assert(0);
    } // switch

  // Loop over cells in fault mesh, computing orientation at each vertex
  const ALE::Obj<sieve_type>& sieve = (*_faultMesh)->getSieve();
  assert(!sieve.isNull());

  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = _quadrature->vertices();
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array vertexOrientation(orientationSize);
  double_array cellVertices(numBasis*spaceDim);

  const ALE::Obj<Mesh::label_sequence>& cellsFault = 
    (*_faultMesh)->heightStratum(0);
  assert(!cellsFault.isNull());
  const Mesh::label_sequence::iterator cellsFaultBegin = cellsFault->begin();
  const Mesh::label_sequence::iterator cellsFaultEnd = cellsFault->end();
  for (Mesh::label_sequence::iterator c_iter=cellsFaultBegin;
       c_iter != cellsFaultEnd;
       ++c_iter) {
    mesh->restrict(coordinates, *c_iter, 
		   &cellVertices[0], cellVertices.size());

    // Compute orientation at each vertex
    int iBasis = 0;
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    assert(!cone.isNull());
    const sieve_type::traits::coneSequence::iterator vBegin = cone->begin();
    const sieve_type::traits::coneSequence::iterator vEnd = cone->end();
    for(sieve_type::traits::coneSequence::iterator v_iter=vBegin;
	v_iter != vEnd;
	++v_iter, ++iBasis) {
      // Compute Jacobian and determinant of Jacobian at vertex
      double_array vertex(&verticesRef[iBasis*cellDim], cellDim);
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, vertex);

      // Compute orientation
      orientFn(&vertexOrientation, jacobian, jacobianDet, upDir);
      
      // Update orientation
      orientation->updatePoint(*v_iter, &vertexOrientation[0]);
    } // for
  } // for

  // Assemble orientation information
  //orientation->complete();

  // Loop over vertices, make orientation information unit magnitude
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    (*_faultMesh)->depthStratum(0);
  const Mesh::label_sequence::iterator vertFaultBegin = vertices->begin();
  const Mesh::label_sequence::iterator vertFaultEnd = vertices->end();
  double_array vertexDir(orientationSize);
  for (Mesh::label_sequence::iterator v_iter=vertFaultBegin;
       v_iter != vertFaultEnd;
       ++v_iter) {
    const real_section_type::value_type* vertexOrient = 
      orientation->restrictPoint(*v_iter);
    
    assert(cellDim*spaceDim == orientationSize);
    for (int iDim=0, index=0; iDim < cellDim; ++iDim, index+=cellDim) {
      double mag = 0;
      for (int jDim=0; jDim < spaceDim; ++jDim)
	mag *= vertexOrient[index*cellDim+jDim];
      for (int jDim=0; jDim < cellDim; ++jDim)
	vertexDir[index*cellDim+jDim] = vertexOrient[index*cellDim+jDim] / mag;
    } // for
    orientation->updatePoint(*v_iter, &vertexDir[0]);
  } // for

  _constraintVert.clear();
  // Create set of vertices associated with Lagrange multiplier constraints
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    // Vertices for each cohesive cell are in groups of N.
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

  // Create orientation section over vertices associated with Lagrange
  // multiplier constraints
  _orientation = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_orientation.isNull());
  const std::set<Mesh::point_type>::const_iterator vertCohesiveBegin = 
    _constraintVert.begin();
  const std::set<Mesh::point_type>::const_iterator vertCohesiveEnd = 
    _constraintVert.end();
  for (std::set<Mesh::point_type>::const_iterator v_iter=vertCohesiveBegin;
       v_iter != vertCohesiveEnd;
       ++v_iter)
    _orientation->setFiberDimension(*v_iter, orientationSize);
  mesh->allocate(_orientation);
  
  // Transfer orientation information from fault vertices to vertices
  // associated with Lagrange multiplier constraints
  Mesh::label_sequence::iterator vFault_iter=vertFaultBegin;
  for (std::set<Mesh::point_type>::const_iterator vConstraint_iter=vertCohesiveBegin;
       vConstraint_iter != vertCohesiveEnd;
       ++vConstraint_iter, ++vFault_iter) {
    const real_section_type::value_type* vertexOrient = 
      orientation->restrictPoint(*vFault_iter);
    _orientation->updatePoint(*vConstraint_iter, vertexOrient);
  } // for
  
  _eqsrc->initialize(mesh, *_faultMesh, _constraintVert, cs);
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
				const ALE::Obj<real_section_type>& residual,
				topology::FieldsManager* const fields,
				const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());
  assert(0 != _faultMesh);
  assert(!_faultMesh->isNull());

  // Subtract constraint forces (which are disp at the constraint
  // DOF) to residual; contributions are at DOF of normal vertices (i and j)

  // Get cohesive cells
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  // Get section information
  const ALE::Obj<real_section_type>& disp = fields->getHistoryItem(1);
  assert(!disp.isNull());  
  
  // Allocate vector for cell values (if necessary)
  _initCellVector();

  // Loop over cohesive cells
  const int numConstraintVert = _quadrature->numBasis();
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    _resetCellVector();

    // Get values at vertices (want constraint forces in disp vector)
    const real_section_type::value_type* cellDisp = 
      mesh->restrict(disp, *c_iter);

    // Transfer constraint forces to cell's constribution to residual vector
    for (int i=0; i < numConstraintVert; ++i) {
      const double constraintForce = cellDisp[2*numConstraintVert+i];
      _cellVector[                  i] = -constraintForce;
      _cellVector[numConstraintVert+i] = -constraintForce;
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numConstraintVert*2);
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

    // Update residual (replace, do not add)
    mesh->update(residual, *c_iter, _cellVector);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveKin::integrateJacobian(
				    PetscMat* mat,
				    topology::FieldsManager* const fields,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != mat);
  assert(0 != fields);
  assert(!mesh.isNull());
  assert(0 != _faultMesh);
  assert(!_faultMesh->isNull());

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  // Get cohesive cells
  const ALE::Obj<ALE::Mesh::label_sequence>& cellsCohesive = 
    mesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const Mesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const Mesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  // Get section information
  const ALE::Obj<real_section_type>& disp = fields->getHistoryItem(1);
  assert(!disp.isNull());  

  const int cellDim = _quadrature->cellDim();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = cellDim*spaceDim;

  // Allocate matrix for cell values (if necessary)
  _initCellMatrix();

  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3*numConstraintVert; // cohesive cell
  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    _resetCellMatrix();

    // Get orientations at cells vertices (only valid at constraint vertices)
    const real_section_type::value_type* cellOrientation = 
      mesh->restrict(_orientation, *c_iter);

    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int iBasis = 3*iConstraint;
      const int jBasis = 3*iConstraint + 1;
      const int kBasis = 3*iConstraint + 2;

      // Get orientations at constraint vertex (3rd vertex)
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[kBasis*orientationSize];

      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = iBasis*spaceDim+iDim;
	  const int col = kBasis*spaceDim+kDim;
	  _cellMatrix[row*numCorners*spaceDim+col] =
	    constraintOrient[iDim*spaceDim+kDim];
	  _cellMatrix[col*numCorners*spaceDim+row] =
	    _cellMatrix[row*numCorners*spaceDim+col]; // symmetric
	} // for

      // Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = jBasis*spaceDim+jDim;
	  const int col = kBasis*spaceDim+kDim;
	  _cellMatrix[row*numCorners*spaceDim+col] =
	    -constraintOrient[jDim*spaceDim+kDim];
	  _cellMatrix[col*numCorners*spaceDim+row] =
	    _cellMatrix[row*numCorners*spaceDim+col]; // symmetric
	} // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numConstraintVert*spaceDim*spaceDim*4);
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", disp);
    // Update values (do not add)
    err = updateOperator(*mat, mesh, disp, globalOrder,
			 *c_iter, _cellMatrix, INSERT_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Set number of degrees of freedom that are constrained at points
void
pylith::faults::FaultCohesiveKin::setConstraintSizes(
				    const ALE::Obj<real_section_type>& field,
				    const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraintSizes
  /* No DOF are eliminated from the system of equations with the 
   * Lagrange multiplier formulation
   */
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::faults::FaultCohesiveKin::setConstraints(
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraints
  /* No DOF are eliminated from the system of equations with the 
   * Lagrange multiplier formulation
   */
} // setConstraints

// ----------------------------------------------------------------------
// Set field.
void
pylith::faults::FaultCohesiveKin::setField(
				     const double t,
				     const ALE::Obj<real_section_type>& disp,
				     const ALE::Obj<Mesh>& mesh)
{ // setField
  typedef std::set<Mesh::point_type>::const_iterator vert_iterator;

  assert(0 != _eqsrc);

  const ALE::Obj<real_section_type>& slip = _eqsrc->slip(t, _constraintVert);
  assert(!slip.isNull());
  const vert_iterator vBegin = _constraintVert.begin();
  const vert_iterator vEnd = _constraintVert.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter)
    disp->updatePoint(*v_iter, slip->restrictPoint(*v_iter));
} // setField


// End of file 
