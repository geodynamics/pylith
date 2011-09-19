// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveTract.hh" // implementation of object methods
#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::DomainSieveMesh SieveMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::RealSection RealSection;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveSubMesh::order_type,PylithInt> IndicesVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveTract::FaultCohesiveTract(void) : 
  _dbInitial(0)
{ // constructor
  _useLagrangeConstraints = false;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveTract::~FaultCohesiveTract(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::FaultCohesiveTract::deallocate(void)
{ // deallocate
  FaultCohesive::deallocate();

  _dbInitial = 0; // :TODO: Use shared pointer
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void pylith::faults::FaultCohesiveTract::dbInitial(spatialdata::spatialdb::SpatialDB* dbs)
{ // dbInitial
  _dbInitial = dbs;
} // dbInitial

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveTract::initialize(const topology::Mesh& mesh,
					       const PylithScalar upDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != _quadrature);

  delete _faultMesh; _faultMesh = new topology::SubMesh();
  CohesiveTopology::createFaultParallel(_faultMesh, mesh, id(),
    useLagrangeConstraints());

  // Reset fields.
  delete _fields; 
  _fields = 
    new topology::Fields<topology::Field<topology::SubMesh> >(*_faultMesh);

  // Initialize quadrature geometry.
  _quadrature->initializeGeometry();

  // Compute orientation at quadrature points in fault mesh.
  _calcOrientation(upDir);

  // Get initial tractions using a spatial database.
  _getInitialTractions();
  
  // Setup fault constitutive model.
  _initConstitutiveModel();
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveTract::integrateResidual(
			   const topology::Field<topology::Mesh>& residual,
			   const PylithScalar t,
			   topology::SolutionFields* const fields)
{ // integrateResidual
  throw std::logic_error("FaultCohesiveTract::integrateResidual() not implemented.");
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveTract::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const PylithScalar t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  throw std::logic_error("FaultCohesiveTract::integrateJacobian() not implemented.");

  _needNewJacobian = false;
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveTract::verifyConfiguration(
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
pylith::faults::FaultCohesiveTract::vertexField(
				       const char* name,
				       const topology::SolutionFields* fields)
{ // vertexField
  throw std::logic_error("FaultCohesiveTract::vertexField() not implemented.");
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveTract::cellField(
				      const char* name,
				      const topology::SolutionFields* fields)
{ // cellField
  throw std::logic_error("FaultCohesiveTract::cellField() not implemented.");
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveTract::_calcOrientation(const PylithScalar upDir[3])
{ // _calcOrientation
  assert(0 != _fields);
  assert(0 != _quadrature);

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

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
  scalar_array quadPtRef(cellDim);
  const scalar_array& quadPtsRef = _quadrature->quadPtsRef();
  
  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int fiberDim = numQuadPts * orientationSize;
  const int jacobianSize = spaceDim * cellDim;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array orientationQuadPt(orientationSize);
  scalar_array orientationCell(fiberDim);

  // Get sections.
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

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
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(PylithScalar));
      cellGeometry.jacobian(&jacobian, &jacobianDet,
			    coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientationQuadPt, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientationQuadPt /= jacobianDet;

      memcpy(&orientationCell[iQuad*orientationSize], 
	     &orientationQuadPt[0], orientationSize*sizeof(PylithScalar));
    } // for

    orientationSection->updatePoint(*c_iter, &orientationCell[0]);
  } // for

  // debugging
  orientation.view("FAULT ORIENTATION");
} // _calcOrientation

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveTract::_getInitialTractions(void)
{ // _getInitialTractions
  assert(0 != _normalizer);
  assert(0 != _quadrature);

  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar lengthScale = _normalizer->lengthScale();

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  if (0 != _dbInitial) { // Setup initial values, if provided.
    // Create section to hold initial tractions.
    _fields->add("initial traction", "initial_traction");
    topology::Field<topology::SubMesh>& traction = _fields->get("initial traction");
    traction.newSection(topology::FieldBase::CELLS_FIELD, numQuadPts*spaceDim);
    traction.allocate();
    traction.scale(pressureScale);
    traction.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
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
    scalar_array tractionCell(numQuadPts*spaceDim);
    scalar_array quadPtsGlobal(numQuadPts*spaceDim);
    
    // Get sections.
    scalar_array coordinatesCell(numBasis*spaceDim);
    const ALE::Obj<RealSection>& coordinates =
      faultSieveMesh->getRealSection("coordinates");
    assert(!coordinates.isNull());
    RestrictVisitor coordsVisitor(*coordinates, 
				  coordinatesCell.size(), &coordinatesCell[0]);

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
      
      const scalar_array& quadPtsNondim = _quadrature->quadPts();
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
pylith::faults::FaultCohesiveTract::_initConstitutiveModel(void)
{ // _initConstitutiveModel
  // :TODO: ADD STUFF HERE
} // _initConstitutiveModel


// End of file 
