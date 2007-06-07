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
  assert(0 != _faultMesh);
  assert(0 != _eqsrc);
  assert(!_faultMesh->isNull());
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");

  // Allocate section for orientation at all fault vertices
  ALE::Obj<real_section_type> orientation = 
    new real_section_type((*_faultMesh)->comm(), (*_faultMesh)->debug());
  assert(!orientation.isNull());
  const int cellDim = (*_faultMesh)->getDimension();
  const int spaceDim = cs->spaceDim();
  const int orientationSize = cellDim*spaceDim;
  orientation->setFiberDimension((*_faultMesh)->depthStratum(0), 
				 orientationSize);
  (*_faultMesh)->allocate(orientation);
  
  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Set orientation method
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

  // Loop over cells, computing orientation at each vertex in cell
  const ALE::Obj<sieve_type>& sieve = (*_faultMesh)->getSieve();
  assert(!sieve.isNull());
  const ALE::Obj<Mesh::label_sequence>& cells = 
    (*_faultMesh)->heightStratum(0);
  const Mesh::label_sequence::iterator cBegin = cells->begin();
  const Mesh::label_sequence::iterator cEnd = cells->end();
  double_array cellOrientation(_quadrature->numBasis()*orientationSize);
  const int numVertices = _quadrature->numBasis();
  for (Mesh::label_sequence::iterator c_iter=cBegin;
       c_iter != cEnd;
       ++c_iter) {
    // Compute Jacobian at vertices
    // STUFF GOES HERE
    double_array jacobian;
    double_array jacobianDet;

    // Compute weighted orientation of face at vertices (using geometry info)
    orientFn(&cellOrientation, jacobian, jacobianDet, upDir, numVertices);

    // Update orientation section for vertices in cell
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    assert(!cone.isNull());
    const sieve_type::traits::coneSequence::iterator vBegin = cone->begin();
    const sieve_type::traits::coneSequence::iterator vEnd = cone->end();
    int index = 0;
    for(sieve_type::traits::coneSequence::iterator v_iter=vBegin;
	v_iter != vEnd;
	++v_iter)
      orientation->updatePoint(*v_iter, &cellOrientation[index]);
      index += orientationSize;
  } // for

  // Assemble orientation information
  //orientation->complete();

  // Loop over vertices, make orientation information unit magnitude
  const ALE::Obj<Mesh::label_sequence>& vertices = 
    (*_faultMesh)->depthStratum(0);
  const Mesh::label_sequence::iterator vBegin = vertices->begin();
  const Mesh::label_sequence::iterator vEnd = vertices->end();
  double_array vertexDir(orientationSize);
  for (Mesh::label_sequence::iterator v_iter=vBegin;
       v_iter != vEnd;
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

  // Create set of constraint vertices
  std::set<Mesh::point_type> setVert;
  for (Mesh::label_sequence::iterator c_iter=cBegin;
       c_iter != cEnd;
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
    for(int i=0, numConstraintVert = coneSize/3; 
	i < numConstraintVert; 
	++i, ++v_iter)
      setVert.insert(*v_iter);
  } // for

  // Only store orientation information at constraint vertices
  _orientation = 
    new real_section_type((*_faultMesh)->comm(), (*_faultMesh)->debug());
  assert(!_orientation.isNull());
  const std::set<Mesh::point_type>::const_iterator cvBegin = 
    _constraintVert.begin();
  const std::set<Mesh::point_type>::const_iterator cvEnd = 
    _constraintVert.end();
  for (std::set<Mesh::point_type>::const_iterator v_iter=cvBegin;
       v_iter != cvEnd;
       ++v_iter)
    _orientation->setFiberDimension(*v_iter, orientationSize);
  (*_faultMesh)->allocate(_orientation);
  for (std::set<Mesh::point_type>::const_iterator v_iter=cvBegin;
       v_iter != cvEnd;
       ++v_iter) {
    const real_section_type::value_type* vertexOrient = 
      orientation->restrictPoint(*v_iter);
    _orientation->updatePoint(*v_iter, vertexOrient);
  } // for
  
  _eqsrc->initialize(mesh, *_faultMesh, setVert, cs);
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

  // Get cell information
  const ALE::Obj<Mesh::label_sequence>& cells = 
    (*_faultMesh)->heightStratum(0);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cBegin = cells->begin();
  const Mesh::label_sequence::iterator cEnd = cells->end();

  // Get section information
  const ALE::Obj<real_section_type>& disp = fields->getHistoryItem(1);
  assert(!disp.isNull());  

  // Allocate vector for cell values (if necessary)
  _initCellVector();

  // Loop over cohesive cells
  const int numVertices = _quadrature->numBasis();
  const int numConstraintVert = numVertices / 3;
  assert(numVertices == numConstraintVert * 3);
  for (Mesh::label_sequence::iterator c_iter=cBegin;
       c_iter != cEnd;
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

    // Update residual
    mesh->updateAdd(residual, *c_iter, _cellVector);
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

  // Get cell informatino
  const ALE::Obj<Mesh::label_sequence>& cells = 
    (*_faultMesh)->heightStratum(0);
  const Mesh::label_sequence::iterator cBegin = cells->begin();
  const Mesh::label_sequence::iterator cEnd = cells->end();

  // Get section information
  const ALE::Obj<real_section_type>& disp = fields->getHistoryItem(1);
  assert(!disp.isNull());  

  const int cellDim = _quadrature->cellDim();
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = cellDim*spaceDim;

  // Allocate matrix for cell values (if necessary)
  _initCellMatrix();

  const int numVertices = _quadrature->numBasis();
  for (Mesh::label_sequence::iterator c_iter=cBegin;
       c_iter != cEnd;
       ++c_iter) {
    _resetCellMatrix();

    // Get orientations at cells vertices (only valid at constraint vertices)
    const real_section_type::value_type* cellOrientation = 
      mesh->restrict(_orientation, *c_iter);

    const int numConstraintVert = numVertices / 3;
    assert(numVertices == numConstraintVert * 3);
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Get orientations at constraint vertex
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[iConstraint*orientationSize];

      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int iBasis = 3*iConstraint;
      const int jBasis = 3*iConstraint + 1;
      const int kBasis = 3*iConstraint + 2;

      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = iBasis*spaceDim+iDim;
	  const int col = kBasis*spaceDim+kDim;
	  _cellMatrix[row*numVertices*spaceDim+col] =
	    constraintOrient[iDim*spaceDim+kDim];
	  _cellMatrix[col*numVertices*spaceDim+row] =
	    _cellMatrix[row*numVertices*spaceDim+col]; // symmetric
	} // for

      // Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = jBasis*spaceDim+jDim;
	  const int col = kBasis*spaceDim+kDim;
	  _cellMatrix[row*numVertices*spaceDim+col] =
	    -constraintOrient[jDim*spaceDim+kDim];
	  _cellMatrix[col*numVertices*spaceDim+row] =
	    _cellMatrix[row*numVertices*spaceDim+col]; // symmetric
	} // for
    } // for
    PetscErrorCode err = 
      PetscLogFlops(numConstraintVert*spaceDim*spaceDim*4);
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", disp);
    err = updateOperator(*mat, mesh, disp, globalOrder,
			 *c_iter, _cellMatrix, ADD_VALUES);
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
  throw std::logic_error("FaultCohesiveKin::setConstraintSizes() not implemented.");
} // setConstraintSizes

// ----------------------------------------------------------------------
// Set which degrees of freedom are constrained at points in field.
void
pylith::faults::FaultCohesiveKin::setConstraints(
				     const ALE::Obj<real_section_type>& field,
				     const ALE::Obj<ALE::Mesh>& mesh)
{ // setConstraints
  throw std::logic_error("FaultCohesiveKin::setConstraints() not implemented.");
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
