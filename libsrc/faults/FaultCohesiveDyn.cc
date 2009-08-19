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
  } // if
} // _getInitialTractions

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveDyn::_initConstitutiveModel(void)
{ // _initConstitutiveModel
  // :TODO: ADD STUFF HERE
} // _initConstitutiveModel


// End of file 
