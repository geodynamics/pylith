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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _fields(0),
  _bufferVectorField(0),
  _bufferScalarField(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  delete _fields; _fields = 0;
  delete _bufferVectorField; _bufferVectorField = 0;
  delete _bufferScalarField; _bufferScalarField = 0;
  // :TODO: Use shared pointers for earthquake sources
} // destructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrcs(const char* const* names,
					 const int numNames,
					 EqKinSrc** sources,
					 const int numSources)
{ // eqsrcs
  assert(numNames == numSources);

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
					     const double upDir[3],
					     const double normalDir[3],
					     spatialdata::spatialdb::SpatialDB* matDB)
{ // initialize
  assert(0 != upDir);
  assert(0 != normalDir);
  assert(0 != _quadrature);

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  
  delete _faultMesh; _faultMesh = new topology::SubMesh();
  CohesiveTopology::createFaultParallel(_faultMesh, &_cohesiveToFault, 
					mesh, id(), _useLagrangeConstraints());

  delete _fields; 
  _fields = new topology::Fields<topology::Field<topology::SubMesh> >(*_faultMesh);

  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter=_eqSrcs.begin(); 
       s_iter != srcsEnd; 
       ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    src->initialize(*_faultMesh, *_normalizer);
  } // for

  // Allocate slip field
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  _fields->add("slip", "slip");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slip.newSection(vertices, cs->spaceDim());
  slip.allocate();

  // Allocate cumulative slip field
  _fields->add("cumulative slip", "cumulative_slip");
  topology::Field<topology::SubMesh>& cumSlip = _fields->get("cumulative slip");
  cumSlip.newSection(slip);
  cumSlip.allocate();

  // Setup pseudo-stiffness of cohesive cells to improve conditioning
  // of Jacobian matrix
  _calcConditioning(cs, matDB);

  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  _quadrature->computeGeometry(*_faultMesh, cells);

  // Compute orientation at vertices in fault mesh.
  _calcOrientation(upDir, normalDir);

  // Compute tributary area for each vertex in fault mesh.
  _calcArea();

  // Create empty tractions field for change in fault tractions.
  _fields->add("tractions", "tractions_change");
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
  assert(0 != _quadrature);
  assert(0 != _fields);

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

  // Allocate vectors for cell values
  double_array orientationCell(numConstraintVert*orientationSize);
  double_array stiffnessCell(numConstraintVert);
  double_array solutionCell(numCorners*spaceDim);
  double_array residualCell(numCorners*spaceDim);
  double_array areaCell(numConstraintVert);
  double_array areaAssembledCell(numConstraintVert);

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get section information
  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  topology::Mesh::RestrictVisitor orientationVisitor(*orientationSection,
						     orientationCell.size(),
						     &orientationCell[0]);

  const ALE::Obj<RealSection>& stiffnessSection = 
    _fields->get("pseudostiffness").section();
  assert(!stiffnessSection.isNull());
  topology::Mesh::RestrictVisitor stiffnessVisitor(*stiffnessSection,
						   stiffnessCell.size(),
						   &stiffnessCell[0]);

  const ALE::Obj<RealSection>& areaSection = 
    _fields->get("area").section();
  assert(!areaSection.isNull());
  topology::Mesh::RestrictVisitor areaVisitor(*areaSection,
					      areaAssembledCell.size(),
					      &areaAssembledCell[0]);

  topology::Field<topology::Mesh>& solution = fields->solution();
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());  
  topology::Mesh::RestrictVisitor solutionVisitor(*solutionSection,
						  solutionCell.size(),
						  &solutionCell[0]);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						   &residualCell[0]);

  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    areaCell = 0.0;
    residualCell = 0.0;

    // Compute contributory area for cell (to weight contributions)
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt*basis[iQuad*numBasis+iBasis];
	areaCell[iBasis] += dArea;
      } // for
    } // for
        
    // Get orientations at fault cell's vertices.
    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);
    
    // Get pseudo stiffness at fault cell's vertices.
    stiffnessVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, stiffnessVisitor);
    
    // Get area at fault cell's vertices.
    areaVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, areaVisitor);
    
    // Get solution at cohesive cell's vertices.
    solutionVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, solutionVisitor);
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      const double pseudoStiffness = stiffnessCell[iConstraint];
      //assert(areaAssembledCell[iConstraint] > 0);
      const double wt = pseudoStiffness * 
	areaCell[iConstraint] / areaAssembledCell[iConstraint];
      
      // Get orientation at constraint vertex
      const double* orientationVertex = 
	&orientationCell[iConstraint*orientationSize];
      assert(0 != orientationVertex);
      
      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim) {
	for (int kDim=0; kDim < spaceDim; ++kDim)
	  residualCell[indexI*spaceDim+iDim] -=
	    solutionCell[indexK*spaceDim+kDim] * 
	    -orientationVertex[kDim*spaceDim+iDim] * wt;
      } // for
      
	// Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	for (int kDim=0; kDim < spaceDim; ++kDim)
	  residualCell[indexJ*spaceDim+jDim] -=
	    solutionCell[indexK*spaceDim+kDim] * 
	    orientationVertex[kDim*spaceDim+jDim] * wt;
      } // for
    } // for

#if 0 // DEBUGGING
    std::cout << "Updating fault residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < numConstraintVert; ++i) {
      std::cout << "  stif["<<i<<"]: " << stiffnessCell[i] << std::endl;
    }
    for(int i = 0; i < numConstraintVert*spaceDim; ++i) {
      std::cout << "  slip["<<i<<"]: " << cellSlip[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  soln["<<i<<"]: " << solutionCell[i] << std::endl;
    }
    for(int i = 0; i < numCorners*spaceDim; ++i) {
      std::cout << "  v["<<i<<"]: " << residualCell[i] << std::endl;
    }
#endif

    residualVisitor.clear();
    sieveMesh->updateAdd(*c_iter, residualVisitor);
  } // for

  // FIX THIS
  PetscLogFlops(cellsCohesiveSize*numConstraintVert*spaceDim*spaceDim*7);
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateResidualAssembled(
			    const topology::Field<topology::Mesh>& residual,
			    const double t,
			    topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  assert(0 != fields);
  assert(0 != _fields);

  // Cohesive cells with normal vertices i and j, and constraint
  // vertex k make 2 contributions to the residual:
  //
  //   * DOF i and j: internal forces in soln field associated with 
  //                  slip
  //   * DOF k: slip values

  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
  slip.zero();
  if (!_useSolnIncr) {
    // Compute slip field at current time step
    const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
    for (srcs_type::iterator s_iter=_eqSrcs.begin(); 
	 s_iter != srcsEnd; 
	 ++s_iter) {
      EqKinSrc* src = s_iter->second;
      assert(0 != src);
      if (t >= src->originTime())
	src->slip(&slip, t);
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
	src->slipIncr(&slip, t-_dt, t);
    } // for
  } // else

  const int spaceDim = _quadrature->spaceDim();

  // Get sections
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& slipSection = slip.section();
  assert(!slipSection.isNull());
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  SieveSubMesh::renumbering_type& renumbering = 
    faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin; 
       v_iter != verticesEnd;
       ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int vertexFault = renumbering[*v_iter];
      const int vertexMesh = *v_iter;
      const double* slipVertex = slipSection->restrictPoint(vertexFault);
      assert(spaceDim == slipSection->getFiberDimension(vertexFault));
      assert(spaceDim == residualSection->getFiberDimension(vertexMesh));
      assert(0 != slipVertex);
      residualSection->updatePoint(vertexMesh, slipVertex);
    } // if
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator that do not
// require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateJacobianAssembled(
				       topology::Jacobian* jacobian,
				       const double t,
				       topology::SolutionFields* const fields)
{ // integrateJacobianAssembled
  assert(0 != jacobian);
  assert(0 != fields);
  assert(0 != _fields);

  typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveMesh::order_type,PetscInt> visitor_type;

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

  PetscErrorCode err = 0;

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();
  const int cellsCohesiveSize = cellsCohesive->size();

  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim*spaceDim;

  const int numConstraintVert = _quadrature->numBasis();
  const int numCorners = 3*numConstraintVert; // cohesive cell
  double_array matrixCell(numCorners*spaceDim * numCorners*spaceDim);
  double_array orientationCell(numConstraintVert*orientationSize);
  double_array stiffnessCell(numConstraintVert);

  // Get section information
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<RealSection>& solutionSection = fields->solution().section();
  assert(!solutionSection.isNull());  
  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  topology::Mesh::RestrictVisitor orientationVisitor(*orientationSection,
						     orientationCell.size(),
						     &orientationCell[0]);

  const ALE::Obj<RealSection>& stiffnessSection = 
    _fields->get("pseudostiffness").section();
  assert(!stiffnessSection.isNull());
  topology::Mesh::RestrictVisitor stiffnessVisitor(*stiffnessSection,
						   stiffnessCell.size(),
						   &stiffnessCell[0]);

#if 0 // DEBUGGING
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

  const PetscMat jacobianMatrix = jacobian->matrix();
  assert(0 != jacobianMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionSection);
  assert(!globalOrder.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  topology::Mesh::IndicesVisitor jacobianVisitor(*solutionSection,
						 *globalOrder,
			   (int) pow(sieveMesh->getSieve()->getMaxConeSize(),
				     sieveMesh->depth())*spaceDim);

  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];

    matrixCell = 0.0;
    // Get orientations at fault cell's vertices.
    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);

    // Get pseudo stiffness at fault cell's vertices.
    stiffnessVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, stiffnessVisitor);
    
    for (int iConstraint=0; iConstraint < numConstraintVert; ++iConstraint) {
      // Blocks in cell matrix associated with normal cohesive
      // vertices i and j and constraint vertex k
      const int indexI = iConstraint;
      const int indexJ = iConstraint +   numConstraintVert;
      const int indexK = iConstraint + 2*numConstraintVert;

      // Get orientation at constraint vertex
      const double* orientationVertex = 
	&orientationCell[iConstraint*orientationSize];
      assert(0 != orientationVertex);

      const double stiffnessVertex = stiffnessCell[iConstraint];

      // Scale orientation information by pseudo-stiffness to bring
      // constraint forces in solution vector to the same order of
      // magnitude as the displacements to prevent ill-conditioning

      // Entries associated with constraint forces applied at node i
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = indexI*spaceDim+iDim;
	  const int col = indexK*spaceDim+kDim;
	  matrixCell[row*numCorners*spaceDim+col] =
	    -orientationVertex[kDim*spaceDim+iDim]*stiffnessVertex;
	  matrixCell[col*numCorners*spaceDim+row] =
	    -orientationVertex[kDim*spaceDim+iDim];
	} // for

      // Entries associated with constraint forces applied at node j
      for (int jDim=0; jDim < spaceDim; ++jDim)
	for (int kDim=0; kDim < spaceDim; ++kDim) {
	  const int row = indexJ*spaceDim+jDim;
	  const int col = indexK*spaceDim+kDim;
	  matrixCell[row*numCorners*spaceDim+col] =
	    orientationVertex[kDim*spaceDim+jDim]*stiffnessVertex;
	  matrixCell[col*numCorners*spaceDim+row] =
	    orientationVertex[kDim*spaceDim+jDim];
	} // for
    } // for

    // Insert cell contribution into PETSc Matrix
    jacobianVisitor.clear();
    err = updateOperator(jacobianMatrix, *sieveMesh->getSieve(),
			 jacobianVisitor, *c_iter, &matrixCell[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
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
  topology::Field<topology::SubMesh>& cumSlip = _fields->get("cumulative slip");
  topology::Field<topology::SubMesh>& slip = _fields->get("slip");
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
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
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
pylith::faults::FaultCohesiveKin::vertexField(
				  const char* name,
				  const topology::SolutionFields* fields)
{ // vertexField
  assert(0 != _faultMesh);
  assert(0 != _quadrature);
  assert(0 != _normalizer);
  assert(0 != _fields);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& cumSlip = 
      _fields->get("cumulative slip");
    return cumSlip;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection =
      orientationSection->getFibration(0);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    assert(0 != _bufferVectorField);
    _bufferVectorField->copy(dirSection);
    _bufferVectorField->label("strike_dir");
    return *_bufferVectorField;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection =
      orientationSection->getFibration(1);
    _allocateBufferVectorField();
    assert(0 != _bufferVectorField);
    _bufferVectorField->copy(dirSection);
    _bufferVectorField->label("dip_dir");
    return *_bufferVectorField;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
    assert(!orientationSection.isNull());
    const int space = 
      (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    const ALE::Obj<RealSection>& dirSection =
      orientationSection->getFibration(space);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    assert(0 != _bufferVectorField);
    _bufferVectorField->copy(dirSection);
    _bufferVectorField->label("normal_dir");
    return *_bufferVectorField;

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen+1);

    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    return s_iter->second->finalSlip();

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    const std::string value = std::string(name).substr(timeStrLen+1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());
    return s_iter->second->slipTime();

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& solution = fields->solution();
    _allocateBufferVectorField();
    _calcTractionsChange(_bufferVectorField, solution);
    _bufferVectorField->label("traction_change");
    return *_bufferVectorField;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name
	<< "' for fault '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  assert(0 != _bufferScalarField);
  return *_bufferScalarField;
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveKin::cellField(
				      const char* name,
				      const topology::SolutionFields* fields)
{ // cellField
  // Should not reach this point if requested field was found
  std::ostringstream msg;
  msg << "Request for unknown cell field '" << name
      << "' for fault '" << label() << ".";
  throw std::runtime_error(msg.str());

  // Return generic section to satisfy member function definition.
  assert(0 != _bufferScalarField);
  return *_bufferScalarField;
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveKin::_calcOrientation(const double upDir[3],
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
  orientation.newSection(topology::FieldBase::VERTICES_FIELD, orientationSize);
  const ALE::Obj<RealSection>& orientationSection = orientation.section();
  assert(!orientationSection.isNull());
  // Create subspaces for along-strike, up-dip, and normal directions
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    orientationSection->addSpace();
  for (int iDim=0; iDim <= cohesiveDim; ++iDim)
    orientationSection->setFiberDimension(vertices, spaceDim, iDim);
  orientation.allocate();
  orientation.zero();
  
  // Get fault cells (1 dimension lower than top-level cells)
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
    // normal with preferred normal is negative, flip up/down dip direction.
    // If the user gives the correct normal direction, we should end
    // up with left-lateral-slip, reverse-slip, and fault-opening for
    // positive slip values.
    
    assert(vertices->size() > 0);
    orientationSection->restrictPoint(*vertices->begin(), &orientationVertex[0],
				      orientationVertex.size());
				      
    assert(3 == spaceDim);
    double_array normalDirVertex(&orientationVertex[6], 3);
    const double dot = 
      normalDir[0]*normalDirVertex[0] +
      normalDir[1]*normalDirVertex[1] +
      normalDir[2]*normalDirVertex[2];
    if (dot < 0.0)
      for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
	   v_iter != verticesEnd;
	   ++v_iter) {
	orientationSection->restrictPoint(*v_iter, &orientationVertex[0],
					  orientationVertex.size());
	assert(9 == orientationSection->getFiberDimension(*v_iter));
	// Flip up-dip direction
	for (int iDim=3; iDim < 6; ++iDim)
	  orientationVertex[iDim] = -orientationVertex[iDim];
	
	// Update direction
	orientationSection->updatePoint(*v_iter, &orientationVertex[0]);
      } // for

    PetscLogFlops(5 + count * 3);
  } // if

  //orientation.view("ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveKin::_calcConditioning(
				 const spatialdata::geocoords::CoordSys* cs,
				 spatialdata::spatialdb::SpatialDB* matDB)
{ // _calcConditioning
  assert(0 != cs);
  assert(0 != matDB);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  const int spaceDim = cs->spaceDim();

  // Get vertices in fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  // Allocate stiffness field.
  _fields->add("pseudostiffness", "pseudostiffness");
  topology::Field<topology::SubMesh>& stiffness = _fields->get("pseudostiffness");
  stiffness.newSection(topology::FieldBase::VERTICES_FIELD, 1);
  stiffness.allocate();
  stiffness.zero();
  const ALE::Obj<RealSection>& stiffnessSection = stiffness.section();
  assert(!stiffnessSection.isNull());

  // Setup queries of physical properties.
  matDB->open();
  const char* stiffnessVals[] = { "density", "vs" };
  const int numStiffnessVals = 2;
  matDB->queryVals(stiffnessVals, numStiffnessVals);
  
  // Get section containing coordinates of vertices
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Get dimensional scales.
  assert(0 != _normalizer);
  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->lengthScale();

  double_array matprops(numStiffnessVals);
  double_array coordsVertex(spaceDim);
  int count = 0;
  
  // Set values in orientation section.
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter, ++count) {
    coordinates->restrictPoint(*v_iter, &coordsVertex[0], coordsVertex.size());
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);
    int err = matDB->query(&matprops[0], numStiffnessVals, &coordsVertex[0], 
			   coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find material properties at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << matDB->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    
    const double density = matprops[0];
    const double vs = matprops[1];
    const double mu = density * vs*vs;
    const double muN = _normalizer->nondimensionalize(mu, pressureScale);
    stiffnessSection->updatePoint(*v_iter, &muN);
  } // for
  PetscLogFlops(count * 2);

  matDB->close();
} // _calcConditioning

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveKin::_calcArea(void)
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
  double_array verticesCell(numBasis*spaceDim);

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
  area.newSection(topology::FieldBase::VERTICES_FIELD, 1);
  area.allocate();
  area.zero();
  const ALE::Obj<RealSection>& areaSection = area.section();
  assert(!areaSection.isNull());
  topology::Mesh::UpdateAddVisitor areaVisitor(*areaSection, &areaCell[0]);  
  
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Loop over cells in fault mesh, compute area
  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    _quadrature->retrieveGeometry(*c_iter);
    areaCell = 0.0;
    
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

// ----------------------------------------------------------------------
// Compute change in tractions on fault surface using solution.
// NOTE: We must convert vertex labels to fault vertex labels
void
pylith::faults::FaultCohesiveKin::_calcTractionsChange(
			     topology::Field<topology::SubMesh>* tractions,
			     const topology::Field<topology::Mesh>& solution)
{ // _calcTractionsChange
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // Get vertices from mesh of domain.
  const ALE::Obj<SieveMesh>& sieveMesh = solution.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();

  // Get fault vertices
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& fvertices = 
    faultSieveMesh->depthStratum(0);
  const SieveSubMesh::label_sequence::iterator fverticesBegin = fvertices->begin();
  const SieveSubMesh::label_sequence::iterator fverticesEnd = fvertices->end();

  // Get sections.
  const ALE::Obj<RealSection>& stiffnessSection = 
    _fields->get("pseudostiffness").section();
  assert(!stiffnessSection.isNull());
  const ALE::Obj<RealSection>& areaSection = 
    _fields->get("area").section();
  assert(!areaSection.isNull());
  const ALE::Obj<RealSection>& solutionSection = solution.section();
  assert(!solutionSection.isNull());  

  const int numFaultVertices = fvertices->size();
  Mesh::renumbering_type& renumbering = faultSieveMesh->getRenumbering();
  const SieveSubMesh::renumbering_type::const_iterator renumberingEnd =
    renumbering.end();

#if 0 // DEBUGGING, MOVE TO SEPARATE CHECK METHOD
  // Check fault mesh and volume mesh coordinates
  const ALE::Obj<RealSection>& coordinates  = mesh->getRealSection("coordinates");
  const ALE::Obj<RealSection>& fCoordinates = _faultMesh->getRealSection("coordinates");

  for (Mesh::label_sequence::iterator v_iter = vertices->begin(); v_iter != verticesEnd; ++v_iter) {
    if (renumbering.find(*v_iter) != renumberingEnd) {
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
  double_array tractionsVertex(fiberDim);

  // Allocate buffer for tractions field (if nec.).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    tractions->newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    tractions->allocate();
  } // if
  assert(!tractionsSection.isNull());
  tractions->zero();
  
  for (SieveMesh::label_sequence::iterator v_iter = verticesBegin; 
       v_iter != verticesEnd;
       ++v_iter)
    if (renumbering.find(*v_iter) != renumberingEnd) {
      const int vertexMesh = *v_iter;
      const int vertexFault = renumbering[*v_iter];
      assert(fiberDim == solutionSection->getFiberDimension(vertexMesh));
      assert(fiberDim == tractionsSection->getFiberDimension(vertexFault));
      assert(1 == stiffnessSection->getFiberDimension(vertexFault));
      assert(1 == areaSection->getFiberDimension(vertexFault));

      const double* solutionVertex = solutionSection->restrictPoint(vertexMesh);
      assert(0 != solutionVertex);
      const double* stiffnessVertex = stiffnessSection->restrictPoint(vertexFault);
      assert(0 != stiffnessVertex);
      const double* areaVertex = areaSection->restrictPoint(vertexFault);
      assert(0 != areaVertex);

      const double scale = stiffnessVertex[0] / areaVertex[0];
      for (int i=0; i < fiberDim; ++i)
	tractionsVertex[i] = solutionVertex[i] * scale;

      tractionsSection->updatePoint(vertexFault, &tractionsVertex[0]);
    } // if

  PetscLogFlops(numFaultVertices * (1 + fiberDim) );

#if 0 // DEBUGGING
  _faultMesh->view("FAULT MESH");
  _area->view("AREA");
  _pseudoStiffness->view("CONDITIONING");
  (*tractions)->view("TRACTIONS");
#endif
} // _calcTractionsChange

// ----------------------------------------------------------------------
// Allocate buffer for vector field.
void
pylith::faults::FaultCohesiveKin::_allocateBufferVectorField(void)
{ // _allocateBufferVectorField
  if (0 != _bufferVectorField)
    return;

  // Create vector field; use same shape/chart as cumulative slip field.
  assert(0 != _faultMesh);
  assert(0 != _fields);
  _bufferVectorField = new topology::Field<topology::SubMesh>(*_faultMesh);
  const topology::Field<topology::SubMesh>& slip = 
    _fields->get("cumulative slip");
  _bufferVectorField->newSection(slip);
  _bufferVectorField->allocate();
  _bufferVectorField->zero();
} // _allocateBufferVectorField

// ----------------------------------------------------------------------
// Allocate buffer for scalar field.
void
pylith::faults::FaultCohesiveKin::_allocateBufferScalarField(void)
{ // _allocateBufferScalarField
  if (0 != _bufferScalarField)
    return;

  // Create vector field; use same shape/chart as area field.
  assert(0 != _faultMesh);
  assert(0 != _fields);
  _bufferScalarField = new topology::Field<topology::SubMesh>(*_faultMesh);
  const topology::Field<topology::SubMesh>& area = _fields->get("area");
  _bufferScalarField->newSection(area);
  _bufferScalarField->allocate();
  _bufferScalarField->zero();
} // _allocateBufferScalarField


// End of file 
