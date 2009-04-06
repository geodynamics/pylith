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
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/utils/array.hh" // USES double_array
#include <petscmat.h> // USES PETSc Mat

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <Completion.hh> // USES completeSection
#include <Selection.hh> // Algorithms for submeshes

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _fields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  delete _fields; _fields = 0;
  // :TODO: Use shared pointers for earthquake sources
} // destructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrcs(const char** names,
					 EqKinSrc** sources,
					 const int numSources)
{ // eqsrcs
  // :TODO: Use shared pointers for earthquake sources
  _eqSrcs.clear();
  for (int i=0; i < numSources; ++i) {
    if (0 == sources[i])
      throw std::runtime_error("Null earthquake source.");
    _eqSrcs[std::string(names[i])] = sources[i];
  } // for
} // eqsrcs

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const topology::Mesh& mesh,
					     const double_array& upDir,
					     const double_array& normalDir,
					     spatialdata::spatialdb::SpatialDB* matDB)
{ // initialize
  assert(0 != _quadrature);

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");
  if (3 != normalDir.size())
    throw std::runtime_error("Normal direction for fault orientation must be "
			     "a vector with 3 components.");

  CohesiveTopology::createParallel(&_faultMesh, &_cohesiveToFault, mesh, id(),
				   _useLagrangeConstraints());
  //_faultMesh->getLabel("height")->view("Fault mesh height");
  //_faultMesh->view("FAULT MESH");

  delete _fields; 
  _fields = new topology::Fields<topology::Field<topology::SubMesh> >;

  // Setup pseudo-stiffness of cohesive cells to improve conditioning
  // of Jacobian matrix
  _calcConditioning(cs, matDB);

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir, normalDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter=_eqSrcs.begin(); 
       s_iter != srcsEnd; 
       ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    src->initialize(_faultMesh, *_normalizer);
  } // for

  // Allocate slip field
  const ALE::Obj<SubMesh::SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  const ALE::Obj<SubMesh::label_sequence>& vertices =
    _faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  _fields->add("slip");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slip.newSection(vertices, cs->spaceDim());
  slip.allocate();

  // Allocate cumulative slip field
  _fields->add("cumulative slip");
  topology::Field<topology::SubMesh>& cumSlip = _fields->get("cumulative slip");
  cumSlip.newSection(slip);
  cumSlip.allocate();
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that
// require assembly across processors.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != fields);

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make 2 contributions to the residual:
  //
  //   * DOF i and j: internal forces in soln field associated with 
  //                  slip
  //   * DOF k: slip values

  if (!_useSolnIncr)
    return;

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim*spaceDim;
  const int numBasis = _quadrature->numBasis();
  const int numConstraintVert = numBasis;
  const int numCorners = 3*numConstraintVert; // cohesive cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const double_array& basis = _quadrature->basis();
  const double_array& jacobianDet = _quadrature->jacobianDet();

  double_array cellOrientation(numConstraintVert*orientationSize);
  double_array cellSoln(numCorners*spaceDim);
  double_array cellStiffness(numConstraintVert);
  double_array cellResidual(numCorners*spaceDim);
  double_array cellArea(numConstraintVert);
  double_array cellAreaAssembled(numConstraintVert);

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = residual.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get section information
  topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<Mesh::RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());  

  for (Mesh::SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const Mesh::SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    cellArea = 0.0;
    cellResidual = 0.0;

    // Compute contributory area for cell (to weight contributions)
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt*basis[iQuad*numBasis+iBasis];
	cellArea[iBasis] += dArea;
      } // for
    } // for
        
    // Get orientations at fault cell's vertices.
    _faultMesh->restrictClosure(_orientation, c_fault, &cellOrientation[0], 
				cellOrientation.size());
    
    // Get pseudo stiffness at fault cell's vertices.
    _faultMesh->restrictClosure(_pseudoStiffness, c_fault, &cellStiffness[0], 
				cellStiffness.size());
    
    // Get pseudo stiffness at fault cell's vertices.
    _faultMesh->restrictClosure(_area, c_fault, &cellAreaAssembled[0], 
				cellAreaAssembled.size());
    
    // Get solution at cohesive cell's vertices.
    mesh->restrictClosure(solution, *c_iter, &cellSoln[0], cellSoln.size());
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      const double pseudoStiffness = cellStiffness[iConstraint];
      const double wt = pseudoStiffness * 
	cellArea[iConstraint] / cellAreaAssembled[iConstraint];
      
      // Get orientation at constraint vertex
      const real_section_type::value_type* constraintOrient = 
	&cellOrientation[iConstraint*orientationSize];
      assert(0 != constraintOrient);
      
      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim) {
	for (int kDim=0; kDim < spaceDim; ++kDim)
	  cellResidual[indexI*spaceDim+iDim] -=
	    cellSoln[indexK*spaceDim+kDim] * 
	    -constraintOrient[kDim*spaceDim+iDim] * wt;
      } // for
      
	// Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	for (int kDim=0; kDim < spaceDim; ++kDim)
	  cellResidual[indexJ*spaceDim+jDim] -=
	    cellSoln[indexK*spaceDim+kDim] * 
	    constraintOrient[kDim*spaceDim+jDim] * wt;
      } // for
    } // for

#if 0 // DEBUGGING
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

    mesh->updateAdd(residual, *c_iter, &cellResidual[0]);
  } // for

  // FIX THIS
  PetscLogFlops(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*7);
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateResidualAssembled(
				const ALE::Obj<real_section_type>& residual,
				const double t,
				topology::FieldsManager* const fields,
				const ALE::Obj<Mesh>& mesh,
				const spatialdata::geocoords::CoordSys* cs)
{ // integrateResidualAssembled
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make 2 contributions to the residual:
  //
  //   * DOF i and j: internal forces in soln field associated with 
  //                  slip
  //   * DOF k: slip values

  assert(!_slip.isNull());
  _slip->zero();
  if (!_useSolnIncr) {
    // Compute slip field at current time step
    const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
    for (srcs_type::iterator s_iter=_eqSrcs.begin(); 
	 s_iter != srcsEnd; 
	 ++s_iter) {
      EqKinSrc* src = s_iter->second;
      assert(0 != src);
      if (t >= src->originTime())
	src->slip(_slip, t, _faultMesh);
    } // for
  } else {
    // Compute increment of slip field at current time step
    const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
    for (srcs_type::iterator s_iter=_eqSrcs.begin(); 
	 s_iter != srcsEnd; 
	 ++s_iter) {
      EqKinSrc* src = s_iter->second;
      assert(0 != src);
      if (t >= src->originTime())
	src->slipIncr(_slip, t-_dt, t, _faultMesh);
    } // for
  } // else

  const int spaceDim = _quadrature->spaceDim();
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  assert(!vertices.isNull());
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();
  Mesh::renumbering_type& renumbering = _faultMesh->getRenumbering();
  for (Mesh::label_sequence::iterator v_iter=vertices->begin(); 
       v_iter != verticesEnd;
       ++v_iter)
    if (renumbering.find(*v_iter) != renumbering.end()) {
      const int vertexFault = renumbering[*v_iter];
      const int vertexMesh = *v_iter;
      const real_section_type::value_type* slip = 
	_slip->restrictPoint(vertexFault);
      assert(spaceDim == _slip->getFiberDimension(vertexFault));
      assert(spaceDim == residual->getFiberDimension(vertexMesh));
      assert(0 != slip);
      residual->updatePoint(vertexMesh, slip);
    } // if
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateJacobianAssembled(
				    PetscMat* mat,
				    const double t,
				    topology::FieldsManager* const fields,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobianAssembled
  assert(0 != mat);
  assert(0 != fields);
  assert(!mesh.isNull());
  typedef ALE::ISieveVisitor::IndicesVisitor<Mesh::real_section_type,Mesh::order_type,PetscInt> visitor_type;

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  PetscErrorCode err = 0;

  // Get cohesive cells
  const ALE::Obj<Mesh::label_sequence>& cellsCohesive = 
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
  double_array cellStiffness(numConstraintVert);

#if 0
  // Check that fault cells match cohesive cells
  ALE::ISieveVisitor::PointRetriever<sieve_type> cV(std::max(1, mesh->getSieve()->getMaxConeSize()));
  ALE::ISieveVisitor::PointRetriever<sieve_type> cV2(std::max(1, _faultMesh->getSieve()->getMaxConeSize()));
  Mesh::renumbering_type& fRenumbering = _faultMesh->getRenumbering();
  const int rank = mesh->commRank();

  for (Mesh::label_sequence::iterator c_iter = cellsCohesiveBegin; c_iter != cellsCohesiveEnd; ++c_iter) {
    mesh->getSieve()->cone(*c_iter, cV);
    const int               coneSize  = cV.getSize();
    const Mesh::point_type *cone      = cV.getPoints();
    const int               faceSize  = coneSize / 3;
    const Mesh::point_type  face      = _cohesiveToFault[*c_iter];
    _faultMesh->getSieve()->cone(face, cV2);
    const int               fConeSize = cV2.getSize();
    const Mesh::point_type *fCone     = cV2.getPoints();

    assert(0 == coneSize % faceSize);
    assert(faceSize == fConeSize);
    // Use last vertices (contraints) for fault mesh
    for(int i = 2*faceSize, j = 0; i < 3*faceSize; ++i, ++j) {
      assert(fRenumbering[cone[i]] == fCone[j]);
    }
    cV.clear();
    cV2.clear();
  }
#endif

  const ALE::Obj<Mesh::order_type>& globalOrder = mesh->getFactory()->getGlobalOrder(mesh, "default", solution);
  assert(!globalOrder.isNull());
  visitor_type iV(*solution, *globalOrder, (int) pow(mesh->getSieve()->getMaxConeSize(), mesh->depth())*spaceDim);

  for (Mesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const Mesh::point_type c_fault = _cohesiveToFault[*c_iter];

    cellMatrix = 0.0;
    // Get orientations at fault cell's vertices.
    _faultMesh->restrictClosure(_orientation, c_fault, &cellOrientation[0], 
				cellOrientation.size());

    // Get pseudo stiffness at fault cell's vertices.
    _faultMesh->restrictClosure(_pseudoStiffness, c_fault, &cellStiffness[0], 
				cellStiffness.size());
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
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

    // Insert cell contribution into PETSc Matrix
    err = updateOperator(*mat, *mesh->getSieve(), iV, *c_iter, &cellMatrix[0], INSERT_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
    iV.clear();
  } // for
  PetscLogFlops(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*4);
  _needNewJacobian = false;
} // integrateJacobianAssembled
  
// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveKin::updateStateVars(const double t,
		       topology::SolutionFields* const fields)
{ // updateStateVars
  assert(0 != fields);
  assert(0 != _fields);

  // Update cumulative slip
  topology::Field<topology::SubMesh>& cumSlip = fields->get("cumulative slip");
  topology::Field<topology::SubMesh>& slip = fields->get("slip");
  if (!_useSolnIncr)
    cumSlip.zero();
  cumSlip += slip;
} // updateStateVars

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveKin::verifyConfiguration(
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

  const int numCorners = _quadrature->refGeometry().numCorners();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  for (Mesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = mesh->getNumCellCorners(*c_iter);
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
  assert(0 != _normalizer);

  const int cohesiveDim = _faultMesh->getDimension();
  const int spaceDim = _quadrature->spaceDim();

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    *fieldType = VECTOR_FIELD;
    assert(!_cumSlip.isNull());
    _allocateBufferVertexVector();
    topology::FieldOps::copyValues(_bufferVertexVector, _cumSlip);
    _bufferTmp = _bufferVertexVector;
    scale = _normalizer->lengthScale();
    fiberDim = spaceDim;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    *fieldType = VECTOR_FIELD;
    _bufferTmp = _orientation->getFibration(0);
    scale = 0.0;
    fiberDim = spaceDim;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    *fieldType = VECTOR_FIELD;
    _bufferTmp = _orientation->getFibration(1);
    scale = 0.0;
    fiberDim = spaceDim;

  } else if (0 == strcasecmp("normal_dir", name)) {
    *fieldType = VECTOR_FIELD;
    const int space = 
      (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    _bufferTmp = _orientation->getFibration(space);
    scale = 0.0;
    fiberDim = spaceDim;

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen+1);

    *fieldType = VECTOR_FIELD;
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    _allocateBufferVertexVector();
    topology::FieldOps::copyValues(_bufferVertexVector, 
				   s_iter->second->finalSlip());
    _bufferTmp = _bufferVertexVector;
    scale = _normalizer->lengthScale();
    fiberDim = spaceDim;

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    *fieldType = SCALAR_FIELD;
    const std::string value = std::string(name).substr(timeStrLen+1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    _allocateBufferVertexScalar();
    topology::FieldOps::copyValues(_bufferVertexScalar, 
				   s_iter->second->slipTime());
    _bufferTmp = _bufferVertexScalar;
    scale = _normalizer->timeScale();
    fiberDim = 1;

  } else if (0 == strcasecmp("traction_change", name)) {
    *fieldType = VECTOR_FIELD;
    const ALE::Obj<real_section_type>& solution = fields->getSolution();
    _calcTractionsChange(&_bufferVertexVector, mesh, solution);
    _bufferTmp = _bufferVertexVector;
    scale = _normalizer->pressureScale();
    fiberDim = spaceDim;
    
  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name
	<< "' for fault '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  if (0 != scale) {
    // dimensionalize values
    double_array values(fiberDim);
    const ALE::Obj<Mesh::label_sequence>& vertices = _faultMesh->depthStratum(0);
    assert(!vertices.isNull());
    const Mesh::label_sequence::iterator verticesEnd = vertices->end();
    for (Mesh::label_sequence::iterator v_iter=vertices->begin(); 
	 v_iter != verticesEnd;
	 ++v_iter) {
      assert(fiberDim == _bufferTmp->getFiberDimension(*v_iter));
      _bufferTmp->restrictPoint(*v_iter, &values[0], values.size());
      _normalizer->dimensionalize(&values[0], values.size(), scale);
      _bufferTmp->updatePointAll(*v_iter, &values[0]);
    } // for
  } // if

  return _bufferTmp;
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const ALE::Obj<pylith::real_section_type>&
pylith::faults::FaultCohesiveKin::cellField(VectorFieldEnum* fieldType,
					    const char* name,
					    const ALE::Obj<Mesh>& mesh,
					    topology::FieldsManager* fields)
{ // cellField
  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown cell field '" << name
      << "' for fault '" << label() << ".";
  throw std::runtime_error(msg.str());

  // Return generic section to satisfy member function definition.
  //return _outputCellVector;
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveKin::_calcOrientation(const double_array& upDir,
						   const double_array& normalDir)
{ // _calcOrientation
  assert(!_faultMesh.isNull());

  // Get vertices in fault mesh
  const ALE::Obj<SubMesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const SubMesh::label_sequence::iterator verticesEnd = vertices->end();

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
  _orientation->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), vertices->end()), *std::max_element(vertices->begin(), vertices->end())+1));
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

  // Loop over cohesive cells, computing orientation weighted by
  // jacobian at constraint vertices
  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = cellGeometry.vertices();
  const int jacobianSize = (cohesiveDim > 0) ? spaceDim * cohesiveDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array vertexOrientation(orientationSize);
  double_array faceVertices(numBasis*spaceDim);
  
  // Get fault cells (1 dimension lower than top-level cells)
  const ALE::Obj<SubMesh::label_sequence>& cells = 
    _faultMesh->heightStratum(0);
  assert(!cells.isNull());
  const SubMesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<sieve_type>& sieve = _faultMesh->getSieve();
  assert(!sieve.isNull());
  const int faultDepth = _faultMesh->depth();  // depth of fault cells
  typedef ALE::SieveAlg<Mesh> SieveAlg;

  ALE::ISieveVisitor::NConeRetriever<sieve_type> ncV(*sieve, (size_t) pow(sieve->getMaxConeSize(), std::max(0, _faultMesh->depth())));

  for (SubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    _faultMesh->restrictClosure(coordinates, *c_iter, 
				&faceVertices[0], faceVertices.size());

    ALE::ISieveTraversal<sieve_type>::orientedClosure(*sieve, *c_iter, ncV);
    const int               coneSize = ncV.getSize();
    const Mesh::point_type *cone     = ncV.getPoints();
    
    for(int v = 0; v < coneSize; ++v) {
      // Compute Jacobian and determinant of Jacobian at vertex
      double_array vertex(&verticesRef[v*cohesiveDim], cohesiveDim);
      cellGeometry.jacobian(&jacobian, &jacobianDet, faceVertices, vertex);

      // Compute orientation
      cellGeometry.orientation(&vertexOrientation, jacobian, jacobianDet, 
			       upDir);
      
      // Update orientation
      _orientation->updateAddPoint(cone[v], &vertexOrientation[0]);
    } // for
    ncV.clear();
  } // for

  // Assemble orientation information
  ALE::Completion::completeSectionAdd(_faultMesh->getSendOverlap(),
				      _faultMesh->getRecvOverlap(),
				      _orientation, _orientation);

  // Loop over vertices, make orientation information unit magnitude
  double_array vertexDir(orientationSize);
  int count = 0;
  for (SubMesh::label_sequence::iterator v_iter=vertices->begin();
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
  PetscLogFlops(count * orientationSize * 4);

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
      for (SubMesh::label_sequence::iterator v_iter=vertices->begin();
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

    PetscLogFlops(5 + count * 3);
  } // if

  //_orientation->view("ORIENTATION");
} // _calcOrientation

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
  const ALE::Obj<SubMesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const SubMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  _pseudoStiffness = new real_section_type(_faultMesh->comm(), 
					   _faultMesh->debug());
  assert(!_pseudoStiffness.isNull());
  _pseudoStiffness->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(), 
									     vertices->end()), 
							   *std::max_element(vertices->begin(), 
									     vertices->end())+1));
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

  assert(0 != _normalizer);
  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->lengthScale();

  double_array matprops(numStiffnessVals);
  double_array vCoords(spaceDim);
  int count = 0;
  
  for (SubMesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    coordinates->restrictPoint(*v_iter, &vCoords[0], vCoords.size());
    _normalizer->dimensionalize(&vCoords[0], vCoords.size(), lengthScale);
    int err = matDB->query(&matprops[0], numStiffnessVals, &vCoords[0], 
			   vCoords.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find material properties at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoords[i];
      msg << ") using spatial database " << matDB->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    
    const double density = matprops[0];
    const double vs = matprops[1];
    const double mu = density * vs*vs;
    const double muN = _normalizer->nondimensionalize(mu, pressureScale);
    _pseudoStiffness->updatePoint(*v_iter, &muN);
  } // for
  PetscLogFlops(count * 2);

  matDB->close();
} // _calcConditioning

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveKin::_calcArea(void)
{ // _calcArea
  assert(!_faultMesh.isNull());

  // Get vertices in fault mesh
  const ALE::Obj<SubMesh::label_sequence>& vertices = 
    _faultMesh->depthStratum(0);
  const SubMesh::label_sequence::iterator verticesEnd = vertices->end();
  const int numVertices = vertices->size();

  _area = new real_section_type(_faultMesh->comm(), 
				_faultMesh->debug());
  assert(!_area.isNull());
  _area->setChart(real_section_type::chart_type(*std::min_element(vertices->begin(),
								  vertices->end()), 
						*std::max_element(vertices->begin(), 
								  vertices->end())+1));
  _area->setFiberDimension(vertices, 1);
  _faultMesh->allocate(_area);
  
  // Get fault cells (1 dimension lower than top-level cells)
  const ALE::Obj<SubMesh::label_sequence>& cells = 
    _faultMesh->heightStratum(0);
  assert(!cells.isNull());
  const SubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    _faultMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Containers for area information
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  double jacobianDet = 0;
  double_array cellArea(numBasis);
  double_array cellVertices(numBasis*spaceDim);

  // Loop over cells in fault mesh, compute area
  for(SubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cellsEnd;
      ++c_iter) {
    _quadrature->computeGeometry(_faultMesh, coordinates, *c_iter);
    cellArea = 0.0;
    
    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute area
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt*basis[iQuad*numBasis+iBasis];
	cellArea[iBasis] += dArea;
      } // for
    } // for
    _faultMesh->updateAdd(_area, *c_iter, &cellArea[0]);

    PetscLogFlops( numQuadPts*(1+numBasis*2) );
  } // for

  // Assemble area information
  ALE::Completion::completeSectionAdd(_faultMesh->getSendOverlap(), _faultMesh->getRecvOverlap(), _area, _area);

#if 0 // DEBUGGING
  _area->view("AREA");
  _faultMesh->getSendOverlap()->view("Send fault overlap");
  _faultMesh->getRecvOverlap()->view("Receive fault overlap");
#endif
} // _calcArea

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
// NOTE: We must convert vertex labels to fault vertex labels
void
pylith::faults::FaultCohesiveKin::_calcTractionsChange(
				 ALE::Obj<real_section_type>* tractions,
				 const ALE::Obj<Mesh>& mesh,
				 const ALE::Obj<real_section_type>& solution)
{ // _calcTractionsChange
  assert(0 != tractions);
  assert(!mesh.isNull());
  assert(!solution.isNull());
  assert(!_faultMesh.isNull());
  assert(!_pseudoStiffness.isNull());
  assert(!_area.isNull());

  const ALE::Obj<Mesh::label_sequence>& vertices = 
    mesh->depthStratum(0);
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const ALE::Obj<Mesh::label_sequence>& fvertices = 
    _faultMesh->depthStratum(0);
  const Mesh::label_sequence::iterator fverticesEnd = fvertices->end();
  const int numFaultVertices = fvertices->size();
  Mesh::renumbering_type& renumbering = _faultMesh->getRenumbering();

#if 0 // MOVE TO SEPARATE CHECK METHOD
  // Check fault mesh and volume mesh coordinates
  const ALE::Obj<real_section_type>& coordinates  = mesh->getRealSection("coordinates");
  const ALE::Obj<real_section_type>& fCoordinates = _faultMesh->getRealSection("coordinates");

  for (Mesh::label_sequence::iterator v_iter = vertices->begin(); v_iter != verticesEnd; ++v_iter) {
    if (renumbering.find(*v_iter) != renumbering.end()) {
      const int     v    = *v_iter;
      const int     dim  = coordinates->getFiberDimension(*v_iter);
      const double *a    = coordinates->restrictPoint(*v_iter);
      const int     fv   = renumbering[*v_iter];
      const int     fDim = fCoordinates->getFiberDimension(fv);
      const double *fa   = fCoordinates->restrictPoint(fv);

      if (dim != fDim) throw ALE::Exception("Coordinate fiber dimensions do not match");
      for(int d = 0; d < dim; ++d) {
        if (a[d] != fa[d]) throw ALE::Exception("Coordinate values do not match");
      }
    }
  }
#endif

  // Fiber dimension of tractions matches spatial dimension.
  const int fiberDim = _quadrature->spaceDim();
  double_array tractionValues(fiberDim);

  // Allocate buffer for tractions field (if nec.).
  if (tractions->isNull() ||
      fiberDim != (*tractions)->getFiberDimension(*fvertices->begin())) {
    *tractions = new real_section_type(_faultMesh->comm(), _faultMesh->debug());
    (*tractions)->setChart(real_section_type::chart_type(
		   *std::min_element(fvertices->begin(), fvertices->end()), 
		   *std::max_element(fvertices->begin(), fvertices->end())+1));
    (*tractions)->setFiberDimension(fvertices, fiberDim);
    _faultMesh->allocate(*tractions);
    assert(!tractions->isNull());
  } // if
  
  for (Mesh::label_sequence::iterator v_iter = vertices->begin(); 
       v_iter != verticesEnd;
       ++v_iter)
    if (renumbering.find(*v_iter) != renumbering.end()) {
      const int vertexMesh = *v_iter;
      const int vertexFault = renumbering[*v_iter];
      assert(fiberDim == solution->getFiberDimension(vertexMesh));
      assert(fiberDim == (*tractions)->getFiberDimension(vertexFault));
      assert(1 == _pseudoStiffness->getFiberDimension(vertexFault));
      assert(1 == _area->getFiberDimension(vertexFault));

      const real_section_type::value_type* solutionValues =
	solution->restrictPoint(vertexMesh);
      assert(0 != solutionValues);
      const real_section_type::value_type* pseudoStiffValue = 
	_pseudoStiffness->restrictPoint(vertexFault);
      assert(0 != _pseudoStiffness);
      const real_section_type::value_type* areaValue = 
	_area->restrictPoint(vertexFault);
      assert(0 != _area);

      const double scale = pseudoStiffValue[0] / areaValue[0];
      for (int i=0; i < fiberDim; ++i)
	tractionValues[i] = solutionValues[i] * scale;

      (*tractions)->updatePoint(vertexFault, &tractionValues[0]);
    } // if

  PetscLogFlops(numFaultVertices * (1 + fiberDim) );

#if 0 // DEBUGGING
  _faultMesh->view("FAULT MESH");
  _area->view("AREA");
  _pseudoStiffness->view("CONDITIONING");
  (*tractions)->view("TRACTIONS");
#endif
} // _calcTractionsChange


// End of file 
