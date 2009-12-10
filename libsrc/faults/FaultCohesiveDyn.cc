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
#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
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
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) : 
  _dbInitial(0)
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

  _dbInitial = 0; // :TODO: Use shared pointer
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void pylith::faults::FaultCohesiveDyn::dbInitial(spatialdata::spatialdb::SpatialDB* dbs)
{ // dbInitial
  _dbInitial = dbs;
} // dbInitial

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

  delete _faultMesh; _faultMesh = new topology::SubMesh();
  CohesiveTopology::createFaultParallel(_faultMesh, &_cohesiveToFault, 
					mesh, id(), useLagrangeConstraints());

  // Reset fields.
  delete _fields; 
  _fields = 
    new topology::Fields<topology::Field<topology::SubMesh> >(*_faultMesh);

  // Initialize quadrature geometry.
  _quadrature->initializeGeometry();

  // Compute orientation at quadrature points in fault mesh.
  _calcOrientation(upDir, normalDir);

  // Get initial tractions using a spatial database.
  _getInitialTractions();
  
  // Setup fault constitutive model.
  _initConstitutiveModel();
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
  assert(0 != _faultMesh);
  assert(0 != _fields);

  // debugging
  residual.view("RESIDUAL BEFORE FRICTION");

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numCornersCohesive = 2*numBasis;
  const int orientationSize = spaceDim*spaceDim;

  // Get cohesive cells.
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cellsCohesive->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cellsCohesive->end();

  // Get fault mesh.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get sections.

  // Get residual.
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  double_array residualCell(numCornersCohesive*spaceDim);
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						   &residualCell[0]);
  double_array forcesCurrentCell(numCornersCohesive*spaceDim);
  RestrictVisitor forcesCurrentVisitor(*residualSection,
				       forcesCurrentCell.size(), &forcesCurrentCell[0]);
  
  // Get displacements.
  double_array dispCell(numCornersCohesive*spaceDim);
  const ALE::Obj<RealSection>& dispSection = fields->get("disp(t)").section();
  assert(!dispSection.isNull());
  RestrictVisitor dispVisitor(*dispSection,
			      dispCell.size(), &dispCell[0]);
  double_array dispIncrCell(numCornersCohesive*spaceDim);
  const ALE::Obj<RealSection>& dispIncrSection = fields->get("dispIncr(t->t+dt)").section();
  assert(!dispIncrSection.isNull());
  RestrictVisitor dispIncrVisitor(*dispIncrSection,
				  dispIncrCell.size(), &dispIncrCell[0]);
  double_array dispTpdtCell(numCornersCohesive*spaceDim);
  
  // Get initial tractions (if exist).
  double_array tractionInitialFault(numQuadPts*spaceDim);
  const ALE::Obj<RealSection>& tractionInitialSection = (0 != _dbInitial) ?
    _fields->get("initial traction").section() : 0;

  // Get orientation.
  double_array orientationCell(numQuadPts*orientationSize);
  const ALE::Obj<RealSection>& orientationSection = _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  // Values defined at quadrature points.
  double_array tractionsCurrentCell(numQuadPts*spaceDim);
  double_array tractionCell(numQuadPts*spaceDim); // friction
  double_array slip(spaceDim); // slip at quad pt
  double_array tractionsCurrentFault(spaceDim); // current tractions at quad pt
  double_array tractionFault(spaceDim); // friction at quad pt

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
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

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get displacements at vertices of cohesive cell.
    dispVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispVisitor);
    dispIncrVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispIncrVisitor);

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtCell = dispCell + dispIncrCell;

    // Get current forces at vertices of cohesive cell (current residual).
    forcesCurrentVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, forcesCurrentVisitor);

    // Get initial tractions for fault cell.
    if (0 != _dbInitial) {
      assert(!tractionInitialSection.isNull());
      tractionInitialSection->restrictPoint(c_fault, &tractionInitialFault[0],
					    tractionInitialFault.size());
    } // if

    // Get fault orientation at quadrature points.
    orientationSection->restrictPoint(c_fault, &orientationCell[0],
				      orientationCell.size());

    tractionCell = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      // wt is also area associated with quadrature point
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];

      // Compute slip at current quad pt in global coordinate system.
      // In: dispTpdtCell [2*numBasis*spaceDim] (negative side then positive side??)
      //     basis [numQuadpts*numBasis]
      // Out: slipGlobal [spaceDim]
      // Use basis functions to compute displacement at quadrature point and
      // then difference displacements to get slip.
      // ADD STUFF HERE
      double_array dispQuadPt(spaceDim*2);
      double_array slipGlobal(spaceDim);
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
	for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	  dispQuadPt[iSpace] += basis[iQuad*numBasis+iBasis]
	    *dispTpdtCell[iBasis*spaceDim+iSpace];
	  dispQuadPt[spaceDim+iSpace] += basis[iQuad*numBasis+iBasis]
	    *dispTpdtCell[(iBasis+numBasis)*spaceDim+iSpace];
	}
      }
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	slipGlobal[iSpace] = dispQuadPt[spaceDim+iSpace]
	  -dispQuadPt[iSpace];
      }
	 
      // Compute slip in fault orientation.
      // In: slipGlobal [spaceDim]
      //     orientationCell [numQuadPts*spaceDim*spaceDim, iQuadPt*spaceDim*spaceDim+iDim*spaceDim+jDim]
      // Out: slipFault [spaceDim]
      // Use orientation to rotate from global to fault orientation.
      // ADD STUFF HERE
      double_array slipFault(spaceDim);
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	for (int jSpace=0; jSpace < spaceDim; ++jSpace) {
	  slipFault[iSpace] += slipGlobal[jSpace]
	    *orientationCell[iQuad*spaceDim*spaceDim+iSpace*spaceDim+jSpace];
	}
      }
	       

      // Compute traction from current deformation in global coordinate system. (tractionCurrentGlobal)
      // In: forcesCurrentCell [2*numBasis*spaceDim] (negative side then positive side)
      // Out: tractionCurrentGlobal [spaceDim]
      // Use basis functions to compute forces at quadrature point, then difference to get relative forces,
      // and divide by area to get traction vector.
      // ADD STUFF HERE
      double_array forcesQuadPt(2*spaceDim);
      double_array tractionCurrentGlobal(spaceDim);
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
	for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	  forcesQuadPt[iSpace] += basis[iQuad*numBasis+iBasis]
	    *forcesCurrentCell[iBasis*spaceDim+iSpace];
	  forcesQuadPt[spaceDim+iSpace] += basis[iQuad*numBasis+iBasis]
	    *forcesCurrentCell[(iBasis+numBasis)*spaceDim+iSpace];
	}
      }
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	tractionCurrentGlobal[iSpace] = forcesQuadPt[spaceDim+iSpace]
	  -forcesQuadPt[iSpace];
	tractionCurrentGlobal[iSpace] /= wt;
      }
      

      // Compute traction in fault orientation.
      // In: tractionCurrentGlobal [spaceDim]
      // Out: tractionCurrentFault [spaceDim]
      // Use orientation to rotate from global to fault orientation.
      // ADD STUFF HERE
      double_array tractionCurrentFault(spaceDim);
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	for (int jSpace=0; jSpace < spaceDim; ++jSpace) {
	  tractionCurrentFault[iSpace] += tractionCurrentGlobal[jSpace]
	    *orientationCell[iQuad*spaceDim*spaceDim+iSpace*spaceDim+jSpace];
	}
      }

      // Compute total traction (initial + current).
      // In: tractionCurrentFault [spaceDim]
      //     tractionInitialFault [numQuadPts*spaceDim]
      // Out: tractionTotalFault [spaceDim]
      // ADD STUFF HERE
      double_array tractionTotalFault(spaceDim); 
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	tractionTotalFault[iSpace] += tractionCurrentFault[iSpace]
	  +tractionInitialFault[iQuad*spaceDim+iSpace];
      }

      // Compute traction ("friction") using fault constitutive model in fault orientation.
      // In: slipFault [spaceDim]
      //     tractionCurrentFault [spaceDim]
      //     tractionTotalFault [spaceDim]
      // Out: frictionFault [spaceDim]
      // BEGIN TEMPORARY
      // Simple fault constitutive model with static friction.
      const double mu = 0.7;
      // ADD STUFF HERE
      double_array frictionFault(spaceDim);
      frictionFault = 0.0;
    switch (spaceDim)
      { // switch
      case 1 : {
	if (tractionTotalFault[0] < 0) {
	  frictionFault[0] = tractionCurrentFault[0];
	    }
	break;
      } // case 1
      case 2 : {
	if (tractionTotalFault[1] < 0) {
	  frictionFault[1] = tractionCurrentFault[1];
	  frictionFault[0] = -mu * tractionTotalFault[1];
	  if (frictionFault[0] > tractionCurrentFault[0])
	    frictionFault[0] = tractionCurrentFault[0];
	}
	break;
      } // case 2
      case 3 : {
	if (tractionTotalFault[2] < 0) {
	  frictionFault[2] = tractionCurrentFault[2];
	  frictionFault[1] = -mu * tractionTotalFault[2] * tractionTotalFault[1] / sqrt(pow(tractionTotalFault[1],2) +pow(tractionTotalFault[0],2));
	  frictionFault[0] = -mu * tractionTotalFault[2] * tractionTotalFault[0] / sqrt(pow(tractionTotalFault[1],2) +pow(tractionTotalFault[0],2));
	
	  if (frictionFault[0] > tractionCurrentFault[0])
	    frictionFault[0] = tractionCurrentFault[0];
	
	  if (frictionFault[1] > tractionCurrentFault[1])
	    frictionFault[1] = tractionCurrentFault[1];
	}
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Friction.");
      } // switch

      
      // END TEMPORARY

      // If normal traction is negative (compression), prevent
      // interpenetration by setting traction to exactly counteract
      // current forces acting normal to fault.

      // Compute traction associated with "friction" in global coordinate system. (tractionCell)
      // In: frictionFault [spaceDim]
      // Out: tractionCell [numQuadPts*spaceDim]
      // Use orientation to rotate from global to fault orientation.
      // ADD STUFF HERE
      for (int iSpace=0; iSpace < spaceDim; ++iSpace) {
	for (int jSpace=0; jSpace < spaceDim; ++jSpace) {
	  tractionCell[iQuad*spaceDim+iSpace] += frictionFault[jSpace]
	    *orientationCell[iQuad*spaceDim*spaceDim+iSpace*spaceDim+jSpace];
	}
      }
      tractionCell /= 2.0;
      
    std::cout << " wt: " << wt
	      << " dispTpdtCell (-): (" << dispTpdtCell[0*spaceDim+0] << "," << dispTpdtCell[0*spaceDim+1] << ")\n"
	      << "                   (" << dispTpdtCell[1*spaceDim+0] << "," << dispTpdtCell[1*spaceDim+1] << ")\n"
	      << " dispTpdtCell (+): (" << dispTpdtCell[2*spaceDim+0] << "," << dispTpdtCell[2*spaceDim+1] << ")\n"
	      << "                   (" << dispTpdtCell[3*spaceDim+0] << "," << dispTpdtCell[3*spaceDim+1] << ")\n"
	      << " dispQuadPt (-): (" << dispQuadPt[0] << "," << dispQuadPt[1] << ")\n"
	      << " dispQuadPt (+): (" << dispQuadPt[spaceDim+0] << "," << dispQuadPt[spaceDim+1] << ")\n"
	      << " slipGlobal: (" << slipGlobal[0] << "," << slipGlobal[1] << ")\n"
	      << " slipFault:  (" << slipFault[0] << "," << slipFault[1] << ")\n"
	      << " forcesQuadPt (-): (" << forcesQuadPt[0] << "," << forcesQuadPt[1] << ")\n"
	      << " forcesQuadPt (+): (" << forcesQuadPt[spaceDim+0] << "," << forcesQuadPt[spaceDim+1] << ")\n"
	      << " tractionCurrentGlobal: (" << tractionCurrentGlobal[0] << "," << tractionCurrentGlobal[1] << ")\n"
	      << " tractionCurrentFault: (" << tractionCurrentFault[0] << "," << tractionCurrentFault[1] << ")\n"
	      << " tractionTotalFault: (" << tractionTotalFault[0] << "," << tractionTotalFault[1] << ")\n"
	      << " frictionFault: (" << frictionFault[0] << "," << frictionFault[1] << ")\n"
	      << " tractionCell: (" << tractionCell[iQuad*spaceDim+0] << "," << tractionCell[iQuad*spaceDim+1] << ")\n"
	      << std::endl;


      // Compute action for dynamic fault term
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
	    // :TODO: action for each side of the fault
            residualCell[iBasis*spaceDim+iDim] += 
	      tractionCell[iQuad*spaceDim+iDim] * valIJ;
	  residualCell[(iBasis+numBasis)*spaceDim+iDim] += 
	      -tractionCell[iQuad*spaceDim+iDim] * valIJ;
	  } // for
        } // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(0)); // :TODO: Count number of operations
  } // for

  // debugging
  residual.view("RESIDUAL AFTER FRICTION");
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveDyn::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const double t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  _needNewJacobian = false;
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

  const int numCorners = _quadrature->refGeometry().numCorners();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieveMesh->getNumCellCorners(*c_iter);
    if (2*numCorners != cellNumCorners) {
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
pylith::faults::FaultCohesiveDyn::_calcOrientation(const double upDir[3],
						   const double normalDir[3])
{ // _calcOrientation
  assert(0 != _fields);
  assert(0 != _quadrature);

  double_array up(upDir, 3);

  // Get 'fault' cells.
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Quadrature related values.
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  double_array quadPtRef(cellDim);
  const double_array& quadPtsRef = _quadrature->quadPtsRef();
  
  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int fiberDim = numQuadPts * orientationSize;
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientationQuadPt(orientationSize);
  double_array orientationCell(fiberDim);

  // Get sections.
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  // :TODO: Use spaces to create subsections like in FaultCohesiveKin.
  _fields->add("orientation", "orientation", 
	       topology::FieldBase::CELLS_FIELD, fiberDim);
  topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");
  orientation.allocate();
  const ALE::Obj<RealSection>& orientationSection = orientation.section();
  assert(!orientationSection.isNull());

  // Loop over cells in fault mesh and compute orientations.
  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    // Compute geometry information for current cell
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);

    // Reset orientation to zero.
    orientationCell = 0.0;

    // Compute orientation at each quadrature point of current cell.
    for (int iQuad=0, iRef=0, iSpace=0; 
	 iQuad < numQuadPts;
	 ++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Reset orientation at quad pt to zero.
      orientationQuadPt = 0.0;

      // Compute Jacobian and determinant at quadrature point, then get
      // orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet,
			    coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientationQuadPt, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientationQuadPt /= jacobianDet;

      memcpy(&orientationCell[iQuad*orientationSize], 
	     &orientationQuadPt[0], orientationSize*sizeof(double));
    } // for

    orientationSection->updatePoint(*c_iter, &orientationCell[0]);
  } // for

  // debugging
  orientation.view("FAULT ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveDyn::_getInitialTractions(void)
{ // _getInitialTractions
  assert(0 != _normalizer);
  assert(0 != _quadrature);

  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->lengthScale();

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  if (0 != _dbInitial) { // Setup initial values, if provided.
    // Create section to hold initial tractions.
    _fields->add("initial traction", "initial_traction");
    topology::Field<topology::SubMesh>& traction = _fields->get("initial traction");
    traction.scale(pressureScale);
    traction.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    traction.newSection(topology::FieldBase::CELLS_FIELD, numQuadPts*spaceDim);
    traction.allocate();
    const ALE::Obj<RealSection>& tractionSection = traction.section();
    assert(!tractionSection.isNull());

    _dbInitial->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbInitial->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbInitial->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-leftlateral",
				    "traction-shear-updip",
				    "traction-normal"};
	_dbInitial->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch

    // Get 'fault' cells.
    const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
    assert(!faultSieveMesh.isNull());
    const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
      faultSieveMesh->heightStratum(0);
    assert(!cells.isNull());
    const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
    const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

    const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
    const int numBasis = _quadrature->numBasis();
    const int numQuadPts = _quadrature->numQuadPts();
    const int spaceDim = _quadrature->spaceDim();
  
    // Containers for database query results and quadrature coordinates in
    // reference geometry.
    double_array tractionCell(numQuadPts*spaceDim);
    double_array quadPtsGlobal(numQuadPts*spaceDim);
    
    // Get sections.
    double_array coordinatesCell(numBasis*spaceDim);
    const ALE::Obj<RealSection>& coordinates =
      faultSieveMesh->getRealSection("coordinates");
    assert(!coordinates.isNull());
    topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						  coordinatesCell.size(),
						  &coordinatesCell[0]);

    const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
    
    // Compute quadrature information
    
    // Loop over cells in boundary mesh and perform queries.
    for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
	 c_iter != cellsEnd;
	 ++c_iter) {
      // Compute geometry information for current cell
      coordsVisitor.clear();
      faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
      _quadrature->computeGeometry(coordinatesCell, *c_iter);
      
      const double_array& quadPtsNondim = _quadrature->quadPts();
      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				  lengthScale);
      
      tractionCell = 0.0;
      for (int iQuad=0, iSpace=0; 
	   iQuad < numQuadPts;
	   ++iQuad, iSpace+=spaceDim) {
	const int err = _dbInitial->query(&tractionCell[iQuad*spaceDim], spaceDim,
					  &quadPtsGlobal[iSpace], spaceDim, cs);
	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial tractions at (";
	  for (int i=0; i < spaceDim; ++i)
	    msg << " " << quadPtsGlobal[i+iSpace];
	  msg << ") for dynamic fault interface " << label() << "\n"
	      << "using spatial database " << _dbInitial->label() << ".";
	  throw std::runtime_error(msg.str());
	} // if
	
      } // for
      _normalizer->nondimensionalize(&tractionCell[0], tractionCell.size(),
				     pressureScale);
      
      // Update section
      assert(tractionCell.size() == tractionSection->getFiberDimension(*c_iter));
      tractionSection->updatePoint(*c_iter, &tractionCell[0]);
    } // for
    
    _dbInitial->close();

    // debugging
    traction.view("INITIAL TRACTIONS");
  } // if
} // _getInitialTractions

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveDyn::_initConstitutiveModel(void)
{ // _initConstitutiveModel
  // :TODO: ADD STUFF HERE
} // _initConstitutiveModel


// End of file 
