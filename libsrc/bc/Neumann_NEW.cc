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

#include "Neumann_NEW.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann_NEW::Neumann_NEW(void) :
  _db(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann_NEW::~Neumann_NEW(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann_NEW::deallocate(void)
{ // deallocate
  _db = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann_NEW::initialize(const topology::Mesh& mesh,
				const double upDir[3])
{ // initialize
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  _queryDatabases();

  _paramsLocalToGlobal(upDir);

  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann_NEW::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != _parameters);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  double_array tractionsCell(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  _calculateValue(t);
  const ALE::Obj<SubRealSection>& tractionSection =
    _parameters->get("value").section();
  assert(!tractionSection.isNull());
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::SubMesh::UpdateAddVisitor residualVisitor(*residualSection,
						      &_cellVector[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    subSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    tractionSection->restrictPoint(*c_iter, 
				&tractionsCell[0], tractionsCell.size());

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
            _cellVector[iBasis*spaceDim+iDim] += 
	      tractionsCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for
    // Assemble cell contribution into field
    residualVisitor.clear();
    subSieveMesh->updateAdd(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::Neumann_NEW::integrateJacobian(topology::Jacobian* jacobian,
					   const double t,
					   topology::SolutionFields* const fields)
{ // integrateJacobian
  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann_NEW::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  if (1 == mesh.dimension())
    throw std::runtime_error("Neumann_NEW boundary conditions are not "
			     "implemented for a 1-D mesh.");

  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::Neumann_NEW::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  assert(0 != _parameters);
  assert(0 != name);

  if (0 == strcasecmp(name, "tractions")) {
    return _parameters->get("traction");

    // ADD STUFF HERE

  } else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Neumann_NEW BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return _parameters->get("traction"); // Satisfy method definition
} // cellField

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::Neumann_NEW::_queryDatabases(void)
{ // _queryDatabases
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  const double rateScale = pressureScale / timeScale;

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int fiberDim = numQuadPts * spaceDim;

  delete _parameters; 
  _parameters = 
    new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);

  // Create section to hold time dependent values
  _parameters->add("value", "traction");
  topology::Field<topology::SubMesh>& value = _parameters->get("value");
  value.scale(pressureScale);
  value.vectorFieldType(topology::FieldBase::OTHER);
  value.newSection(topology::FieldBase::CELLS_FIELD, fiberDim);
  value.allocate();

  if (0 != _dbInitial) { // Setup initial values, if provided.
    _parameters->add("initial", "initial_traction");
    topology::Field<topology::SubMesh>& initial = 
      _parameters->get("initial");
    initial.cloneSection(value);
    initial.scale(pressureScale);
    initial.vectorFieldType(topology::FieldBase::MULTI_VECTOR);

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
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbInitial->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann_NEW.");
      } // switch
    _queryDB(&initial, _dbInitial, spaceDim, pressureScale);
    _dbInitial->close();
  } // if

  if (0 != _dbRate) { // Setup rate of change of values, if provided.
    _parameters->add("rate", "traction_rate");
    topology::Field<topology::SubMesh>& rate = 
      _parameters->get("rate");
    rate.cloneSection(value);
    rate.scale(rateScale);
    rate.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    const ALE::Obj<RealSection>& rateSection = rate.section();
    assert(!rateSection.isNull());

    _dbRate->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-rate-normal"};
	_dbRate->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-rate-shear", 
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-rate-shear-horiz",
				    "traction-rate-shear-vert",
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann_NEW.");
      } // switch
    _queryDB(&rate, _dbRate, spaceDim, pressureScale);

    _parameters->add("rate time", "rate_traction_time");
    topology::Field<topology::SubMesh>& rateTime = 
      _parameters->get("rate time");
    rateTime.newSection(rate, 1);
    rateTime.allocate();
    rateTime.scale(timeScale);
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB(&rateTime, _dbRate, 1, timeScale);
    _dbRate->close();
  } // if

  if (0 != _dbChange) { // Setup change of values, if provided.
    _parameters->add("change", "change_traction");
    topology::Field<topology::SubMesh>& change = 
      _parameters->get("change");
    change.cloneSection(value);
    change.scale(pressureScale);
    change.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    const ALE::Obj<RealSection>& changeSection = change.section();
    assert(!changeSection.isNull());

    _dbChange->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbChange->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbChange->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-horiz",
				    "traction-shear-vert",
				    "traction-normal"};
	_dbChange->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in Neumann_NEW.");
      } // switch
    _queryDB(&change, _dbChange, spaceDim, pressureScale);

    _parameters->add("change time", "change_traction_time");
    topology::Field<topology::SubMesh>& changeTime = 
      _parameters->get("change time");
    changeTime.newSection(change, 1);
    changeTime.allocate();
    changeTime.scale(timeScale);
    changeTime.vectorFieldType(topology::FieldBase::SCALAR);

    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB(&changeTime, _dbChange, 1, timeScale);
    _dbChange->close();

    if (0 != _dbTimeHistory)
      _dbTimeHistory->open();
  } // if

} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::Neumann_NEW::_queryDB(topology::Field<topology::SubMesh>* field,
			      spatialdata::spatialdb::SpatialDB* const db,
			      const int querySize,
			      const double scale)
{ // _queryDB
  assert(0 != field);
  assert(0 != db);
  assert(0 != _boundaryMesh);
  assert(0 != _quadrature);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  const int fiberDim = numQuadPts * querySize;
  
  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  double_array valuesCell(querySize);
  double_array quadPtRef(cellDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);
  const double_array& quadPtsRef = _quadrature->quadPtsRef();

  // Get sections.
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    subSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  const ALE::Obj<RealSection>& section = field->section();
  assert(!section.isNull());

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif

  // Loop over cells in boundary mesh and perform queries.
  for (SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
    const double_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    valuesCell = 0.0;
    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts;
	++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      const int err = _db->query(&valuesCell[iQuad*querySize], querySize,
				 &quadPtsGlobal[iSpace], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find values at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << " " << quadPtsGlobal[i+iSpace];
	msg << ") for traction boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      _normalizer->nondimensionalize(&valuesCell[0], valuesCell.size(),
				     scale);

    } // for

    // Update section
    section->updatePoint(*c_iter, &valuesCell[0]);
  } // for
} // _queryDB

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
  pylith::bc::Neumann_NEW::_paramsLocalToGlobal(const double upDir[3])
{ // _paramsLocalToGlobal
  assert(0 != _boundaryMesh);
  assert(0 != _quadrature);

  double_array up(upDir, 3);


  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Quadrature related values.
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cellGeometry.spaceDim();
  const int fiberDim = spaceDim * numQuadPts;
  double_array quadPtRef(cellDim);
  const double_array& quadPtsRef = _quadrature->quadPtsRef();
  
  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const int jacobianSize = spaceDim * cellDim;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);

  // Container for cell tractions rotated to global coordinates.
  double_array initialCellLocal(fiberDim);
  double_array initialCellGlobal(fiberDim);
  double_array rateCellLocal(fiberDim);
  double_array rateCellGlobal(fiberDim);
  double_array changeCellLocal(fiberDim);
  double_array changeCellGlobal(fiberDim);

  // Get sections.
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    subSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordinates, 
						coordinatesCell.size(),
						&coordinatesCell[0]);

  const ALE::Obj<RealSection>& initialSection = (0 != _dbInitial) ?
    _parameters->get("initial").section() : 0;
  const ALE::Obj<RealSection>& rateSection = ( 0 != _dbRate) ?
    _parameters->get("rate").section() : 0;
  const ALE::Obj<RealSection>& changeSection = ( 0 != _dbChange) ?
    _parameters->get("change").section() : 0;

  // Loop over cells in boundary mesh, compute orientations, and then
  // rotate corresponding traction vector from local to global coordinates.
  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif
    // Reset traction vectors
    initialCellLocal = 0.0;
    initialCellGlobal = 0.0;
    rateCellLocal = 0.0;
    rateCellGlobal = 0.0;
    changeCellLocal = 0.0;
    changeCellGlobal = 0.0;

    // Get values for cell from each of the sections
    if (0 != _dbInitial) {
      assert(!initialSection.isNull());
      initialSection->restrictPoint(*c_iter, &initialCellLocal[0],
				    initialCellLocal.size());
    } // if
    if (0 != _dbRate) {
      assert(!rateSection.isNull());
      rateSection->restrictPoint(*c_iter, &rateCellLocal[0], 
				 rateCellLocal.size());
    } // if
    if (0 != _dbChange) {
      assert(!changeSection.isNull());
      changeSection->restrictPoint(*c_iter, &changeCellLocal[0], 
				   changeCellLocal.size());
    } // if

    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts;
	++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Compute Jacobian and determinant at quadrature point, then get
      // orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(double));
#if defined(PRECOMPUTE_GEOMETRY)
      coordsVisitor.clear();
      subSieveMesh->restrictClosure(*c_iter, coordsVisitor);
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet,
			    coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      if (0 != _dbInitial) {
	assert(!initialSection.isNull());
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    initialCellGlobal[iDim+iSpace] +=
	      orientation[jDim*spaceDim+iDim] * initialCellLocal[jDim];
	} // for
      } // if

      if (0 != _dbRate) {
	assert(!rateSection.isNull());
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    rateCellGlobal[iDim+iSpace] +=
	      orientation[jDim*spaceDim+iDim] * rateCellLocal[jDim];
	} // for
      } // if

      if (0 != _dbChange) {
	assert(!changeSection.isNull());
	// Rotate traction vector from local coordinate system to global
	// coordinate system
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  for(int jDim = 0; jDim < spaceDim; ++jDim)
	    changeCellGlobal[iDim+iSpace] +=
	      orientation[jDim*spaceDim+iDim] * changeCellLocal[jDim];
	} // for
      } // if

    } // for
    
    // Update sections
    if (0 != _dbInitial)
      initialSection->updatePoint(*c_iter, &initialCellGlobal[0]);
    if (0 != _dbRate)
      rateSection->updatePoint(*c_iter, &rateCellGlobal[0]);
    if (0 != _dbChange)
      changeSection->updatePoint(*c_iter, &changeCellGlobal[0]);
  } // for
} // paramsLocalToGlobal

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::bc::Neumann_NEW::_calculateValue(const double t)
{ // _calculateValue
  assert(0 != _parameters);
  assert(0 != _boundaryMesh);
  assert(0 != _quadrature);

  const ALE::Obj<RealSection>& valueSection = 
    _parameters->get("value").section();
  assert(!valueSection.isNull());
  valueSection->zero();

  const double timeScale = _getNormalizer().timeScale();

  const ALE::Obj<RealSection>& initialSection = (0 != _dbInitial) ?
    _parameters->get("initial").section() : 0;
  const ALE::Obj<RealSection>& rateSection = ( 0 != _dbRate) ?
    _parameters->get("rate").section() : 0;
  const ALE::Obj<RealSection>& rateTimeSection = (0 != _dbRate) ?
    _parameters->get("rate time").section() : 0;
  const ALE::Obj<RealSection>& changeSection = ( 0 != _dbChange) ?
    _parameters->get("change").section() : 0;
  const ALE::Obj<RealSection>& changeTimeSection = ( 0 != _dbChange) ?
    _parameters->get("change time").section() : 0;

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& subSieveMesh = _boundaryMesh->sieveMesh();
  assert(!subSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    subSieveMesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  const int spaceDim = _quadrature->spaceDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int fiberDim = spaceDim * numQuadPts;

  double_array valuesCell(fiberDim);
  double_array bufferCell(fiberDim);
  for(SieveSubMesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    
    valuesCell = 0.0;
    
    // Contribution from initial value
    if (0 != _dbInitial) {
      assert(!initialSection.isNull());
      initialSection->restrictPoint(*c_iter, 
				    &bufferCell[0], bufferCell.size());
      valuesCell += bufferCell;
    } // if
    
    // Contribution from rate of change of value
    if (0 != _dbRate) {
      assert(!rateSection.isNull());
      assert(!rateTimeSection.isNull());
      double tRate = 0.0;
      
      rateSection->restrictPoint(*c_iter, &bufferCell[0], bufferCell.size());
      rateTimeSection->restrictPoint(*c_iter, &tRate, 1);
      if (t > tRate) { // rate of change integrated over time
	bufferCell *= (t - tRate);
	valuesCell += bufferCell;
      } // if
    } // if
    
    // Contribution from change of value
    if (0 != _dbChange) {
      assert(!changeSection.isNull());
      assert(!changeTimeSection.isNull());
      double tChange = 0.0;

      changeSection->restrictPoint(*c_iter, &bufferCell[0], bufferCell.size());
      changeTimeSection->restrictPoint(*c_iter, &tChange, 1);
      if (t >= tChange) { // change in value over time
	double scale = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = t - tChange;
	  _getNormalizer().dimensionalize(&tDim, 1, timeScale);
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (0 != err) {
	    std::ostringstream msg;
	    msg << "Error querying for time '" << tDim 
		<< "' in time history database "
		<< _dbTimeHistory->label() << ".";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	bufferCell *= scale;
	valuesCell += bufferCell;
      } // if
    } // if

    valueSection->updateAddPoint(*c_iter, &valuesCell[0]);
  } // for
}  // _calculateValue


// End of file 
