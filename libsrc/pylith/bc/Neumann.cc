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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Neumann.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field
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
// Default constructor.
pylith::bc::Neumann::Neumann(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::Neumann::deallocate(void)
{ // deallocate
  BCIntegratorSubMesh::deallocate();
  TimeDependent::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const topology::Mesh& mesh,
				const PylithScalar upDir[3])
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
pylith::bc::Neumann::integrateResidual(const topology::Field<topology::Mesh>& residual,
				       const PylithScalar t,
				       topology::SolutionFields* const fields)
{ // integrateResidual
  assert(_quadrature);
  assert(_boundaryMesh);
  assert(_parameters);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const scalar_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  // Allocate vectors for cell values.
  _initCellVector();
  scalar_array tractionsCell(numQuadPts*spaceDim);

  // Get cell information
  PetscDM subMesh = _boundaryMesh->dmMesh();assert(subMesh);
  PetscIS subpointIS;
  PetscInt cStart, cEnd;
  PetscErrorCode err = 0;
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
  err = DMPlexCreateSubpointIS(subMesh, &subpointIS);CHECK_PETSC_ERROR(err);

  // Get sections
  _calculateValue(t);
  topology::Field<topology::SubMesh>& valueField = _parameters->get("value");
  
  PetscSection residualSection = residual.petscSection();assert(residualSection);
  PetscVec residualVec = residual.localVector();assert(residualVec);
  PetscSection residualSubsection = NULL;
  err = PetscSectionCreateSubmeshSection(residualSection, subpointIS, &residualSubsection);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&subpointIS);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  PetscScalar* valueArray = valueField.getLocalArray();

  // Loop over faces and integrate contribution from each face
  for(PetscInt c = cStart; c < cEnd; ++c) {
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#else
    const PetscScalar *coords;
    PetscInt coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for (PetscInt i=0; i < coordsSize; ++i) {
      coordinatesCell[i] = coords[i];
    } // for
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    const PetscInt voff = valueField.sectionOffset(c);
    assert(numQuadPts*spaceDim == valueField.sectionDof(c));

    // Get cell geometry information that depends on cell
    const scalar_array& basis = _quadrature->basis();
    const scalar_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const PylithScalar valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const PylithScalar valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += valueArray[voff+iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for
    err = DMPlexVecSetClosure(subMesh, residualSubsection, residualVec, c, &_cellVector[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
  valueField.restoreLocalArray(&valueArray);
} // integrateResidual

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  if (1 == mesh.dimension())
    throw std::runtime_error("Neumann boundary conditions are not "
			     "implemented for a 1-D mesh.");

  BCIntegratorSubMesh::verifyConfiguration(mesh);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field for tractions.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::bc::Neumann::cellField(const char* name,
			       topology::SolutionFields* const fields)
{ // cellField
  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "initial_value"))
    return _parameters->get("initial");

  else if (0 == strcasecmp(name, "rate_of_change"))
    return _parameters->get("rate");

  else if (0 == strcasecmp(name, "change_in_value"))
    return _parameters->get("change");

  else if (0 == strcasecmp(name, "rate_start_time"))
    return _parameters->get("rate time");

  else if (0 == strcasecmp(name, "change_start_time"))
    return _parameters->get("change time");

  else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for Neumann BC '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return _parameters->get("traction"); // Satisfy method definition
} // cellField

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::Neumann::_queryDatabases(void)
{ // _queryDatabases
  assert(_quadrature);
  assert(_boundaryMesh);
  
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  const PylithScalar rateScale = pressureScale / timeScale;

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();

  delete _parameters; 
  _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(*_boundaryMesh);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
  topology::Field<topology::SubMesh>& value = _parameters->get("value");
  value.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
  value.scale(pressureScale);
  value.allocate();
  if (_dbInitial) {
    _parameters->add("initial", "initial_traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field<topology::SubMesh>& initial = _parameters->get("initial");
    initial.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    initial.scale(pressureScale);
    initial.allocate();
  }
  if (_dbRate) {
    _parameters->add("rate", "traction_rate", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field<topology::SubMesh>& rate = _parameters->get("rate");
    rate.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    rate.scale(rateScale);
    rate.allocate();
    _parameters->add("rate time", "traction_rate_time", topology::FieldBase::FACES_FIELD, numQuadPts);
    topology::Field<topology::SubMesh>& rateTime = _parameters->get("rate time");
    rateTime.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    rateTime.scale(timeScale);
    rateTime.allocate();
  } // if
  if (_dbChange) {
    _parameters->add("change", "change_traction", topology::FieldBase::FACES_FIELD, numQuadPts*spaceDim);
    topology::Field<topology::SubMesh>& change = _parameters->get("change");
    change.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    change.scale(pressureScale);
    change.allocate();
    _parameters->add("change time", "change_traction_time", topology::FieldBase::FACES_FIELD, numQuadPts);
    topology::Field<topology::SubMesh>& changeTime = _parameters->get("change time");
    changeTime.vectorFieldType(topology::FieldBase::MULTI_VECTOR);
    changeTime.scale(timeScale);
    changeTime.allocate();
  } // if

  if (_dbInitial) { // Setup initial values, if provided.
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
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("initial", _dbInitial, spaceDim, pressureScale);
    _dbInitial->close();
  } // if

  if (_dbRate) { // Setup rate of change of values, if provided.
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
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("rate", _dbRate, spaceDim, rateScale);

    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale);
    _dbRate->close();
  } // if

  if (_dbChange) { // Setup change of values, if provided.
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
	throw std::logic_error("Bad spatial dimension in Neumann.");
      } // switch
    _queryDB("change", _dbChange, spaceDim, pressureScale);

    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale);
    _dbChange->close();

    if (_dbTimeHistory)
      _dbTimeHistory->open();
  } // if

} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::Neumann::_queryDB(const char* name,
			      spatialdata::spatialdb::SpatialDB* const db,
			      const int querySize,
			      const PylithScalar scale)
{ // _queryDB
  assert(name);
  assert(db);
  assert(_boundaryMesh);
  assert(_quadrature);
  assert(_parameters);

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM subMesh = _boundaryMesh->dmMesh();assert(subMesh);
  PetscInt cStart, cEnd;
  PetscErrorCode err = 0;
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  const int cellDim = _quadrature->cellDim() > 0 ? _quadrature->cellDim() : 1;
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = _quadrature->spaceDim();
  
  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valuesCell(numQuadPts*querySize);
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get sections.
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  topology::Field<topology::SubMesh>& valueField = _parameters->get(name);

  const spatialdata::geocoords::CoordSys* cs = _boundaryMesh->coordsys();
  assert(cs);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  // Compute quadrature information
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  _quadrature->computeGeometry(*_boundaryMesh, cells);
#endif
  
  PetscScalar* valueArray = valueField.getLocalArray();

  // Loop over cells in boundary mesh and perform queries.
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#else
    const PetscScalar *coords;
    PetscInt coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for (PetscInt i = 0; i < coordsSize; ++i) {
      coordinatesCell[i] = coords[i];
    } // for
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif
    const scalar_array& quadPtsNondim = _quadrature->quadPts();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);
    
    valuesCell = 0.0;
    for (int iQuad=0, iSpace=0; iQuad < numQuadPts; ++iQuad, iSpace+=spaceDim) {
      const int err = db->query(&valuesCell[iQuad*querySize], querySize,
                                &quadPtsGlobal[iSpace], spaceDim, cs);
      
      if (err) {
        std::ostringstream msg;
        msg << "Could not find values at (";
        for (int i=0; i < spaceDim; ++i)
          msg << " " << quadPtsGlobal[i+iSpace];
        msg << ") for traction boundary condition " << _label << "\n"
            << "using spatial database " << db->label() << ".";
        throw std::runtime_error(msg.str());
      } // if
    } // for
    _normalizer->nondimensionalize(&valuesCell[0], valuesCell.size(), scale);

    // Update section
    const PetscInt voff = valueField.sectionOffset(c);
    const PetscInt vdof = valueField.sectionDof(c);
    assert(numQuadPts*querySize == vdof);
    for(PetscInt d = 0; d < vdof; ++d)
      valueArray[voff+d] = valuesCell[d];
  } // for
  valueField.restoreLocalArray(&valueArray);
} // _queryDB

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
  pylith::bc::Neumann::_paramsLocalToGlobal(const PylithScalar upDir[3])
{ // _paramsLocalToGlobal
  assert(_boundaryMesh);
  assert(_parameters);
  assert(_quadrature);

  scalar_array up(3);
  for (int i=0; i < 3; ++i)
    up[i] = upDir[i];

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM subMesh = _boundaryMesh->dmMesh();assert(subMesh);
  PetscInt cStart, cEnd;
  PetscErrorCode err = 0;
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

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
  const int jacobianSize = spaceDim * cellDim;
  scalar_array jacobian(jacobianSize);
  PylithScalar jacobianDet = 0;
  scalar_array orientation(orientationSize);

  // Get sections.
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  err = DMPlexGetCoordinateSection(subMesh, &coordSection);CHECK_PETSC_ERROR(err);assert(coordSection);
  err = DMGetCoordinatesLocal(subMesh, &coordVec);CHECK_PETSC_ERROR(err);assert(coordVec);

  scalar_array parametersCellLocal(spaceDim);

  topology::Field<topology::SubMesh>& valueField = _parameters->get("value");
  PetscScalar* valueArray = valueField.getLocalArray();

  topology::Field<topology::SubMesh>* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::Field<topology::SubMesh>* rateField = (_dbRate) ? &_parameters->get("rate") : 0;
  topology::Field<topology::SubMesh>* rateTimeField = (_dbRate) ? &_parameters->get("rate time") : 0;
  topology::Field<topology::SubMesh>* changeField = (_dbChange) ? &_parameters->get("change") : 0;
  topology::Field<topology::SubMesh>* changeTimeField = (_dbChange) ? &_parameters->get("change time") : 0;

  PetscScalar* initialArray = (_dbInitial) ? initialField->getLocalArray() : NULL;
  PetscScalar* rateArray = (_dbRate) ? rateField->getLocalArray() : NULL;
  PetscScalar* rateTimeArray = (_dbRate) ? rateTimeField->getLocalArray() : NULL;
  PetscScalar* changeArray = (_dbChange) ? changeField->getLocalArray() : NULL;
  PetscScalar* changeTimeArray = (_dbChange) ? changeTimeField->getLocalArray() : NULL;

  // Loop over cells in boundary mesh, compute orientations, and then
  // rotate corresponding traction vector from local to global coordinates.
  for(PetscInt c = cStart; c < cEnd; ++c) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#else
    const PetscScalar *coords;
    PetscInt coordsSize;
    err = DMPlexVecGetClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for (PetscInt i=0; i < coordsSize; ++i) {
      coordinatesCell[i] = coords[i];
    } // for
    _quadrature->computeGeometry(coordinatesCell, c);
    err = DMPlexVecRestoreClosure(subMesh, coordSection, coordVec, c, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif
    for(int iQuad=0, iRef=0, iSpace=0; iQuad < numQuadPts; ++iQuad, iRef+=cellDim, iSpace+=spaceDim) {
      // Compute Jacobian and determinant at quadrature point, then get orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[iRef], cellDim*sizeof(PylithScalar));
#if defined(PRECOMPUTE_GEOMETRY)
#error("Code for PRECOMPUTE_GEOMETRY not implemented.")
#endif
      cellGeometry.jacobian(&jacobian, &jacobianDet, coordinatesCell, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, up);
      assert(jacobianDet > 0.0);
      orientation /= jacobianDet;

      if (_dbInitial) {
        // Rotate traction vector from local coordinate system to global coordinate system
        PylithScalar *initialLocal = &parametersCellLocal[0];
	assert(initialField);
	const PetscInt ioff = initialField->sectionOffset(c);
	assert(numQuadPts*spaceDim == initialField->sectionDof(c));
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          initialLocal[iDim] = initialArray[ioff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          initialArray[ioff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            initialArray[ioff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * initialLocal[jDim];
	  } // for
        } // for
      } // if

      if (_dbRate) {
        // Rotate traction vector from local coordinate system to global coordinate system
        PylithScalar *rateLocal = &parametersCellLocal[0];
	assert(rateField);
	const PetscInt roff = rateField->sectionOffset(c);
	assert(numQuadPts*spaceDim == rateField->sectionDof(c));
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          rateLocal[iDim] = rateArray[roff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          rateArray[roff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            rateArray[roff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * rateLocal[jDim];
	  } // for
        } // for
      } // if

      if (_dbChange) {
        // Rotate traction vector from local coordinate system to global coordinate system
        PylithScalar *changeLocal = &parametersCellLocal[0];
	assert(changeField);
	const PetscInt coff = changeField->sectionOffset(c);
	assert(numQuadPts*spaceDim == changeField->sectionDof(c));
        for (int iDim = 0; iDim < spaceDim; ++iDim) {
          changeLocal[iDim] = changeArray[coff+iSpace+iDim];
        } // for
        for(int iDim = 0; iDim < spaceDim; ++iDim) {
          changeArray[coff+iSpace+iDim] = 0.0;
          for(int jDim = 0; jDim < spaceDim; ++jDim) {
            changeArray[coff+iSpace+iDim] += orientation[jDim*spaceDim+iDim] * changeLocal[jDim];
	  } // for
        } // for
      } // if
    } // for
  } // for

  valueField.restoreLocalArray(&valueArray);
  if (_dbInitial) {
    assert(initialField);
    initialField->restoreLocalArray(&initialArray);
  } // if
  if (_dbRate) {
    assert(rateField);
    rateField->restoreLocalArray(&rateArray);
    assert(rateTimeField);
    rateTimeField->restoreLocalArray(&rateTimeArray);
  } // if
  if (_dbChange) {
    assert(changeField);
    changeField->restoreLocalArray(&changeArray);
    assert(changeTimeField);
    changeTimeField->restoreLocalArray(&changeTimeArray);
  } // if

} // paramsLocalToGlobal

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::bc::Neumann::_calculateValue(const PylithScalar t)
{ // _calculateValue
  assert(_parameters);
  assert(_boundaryMesh);
  assert(_quadrature);

  const PylithScalar timeScale = _getNormalizer().timeScale();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  PetscDM subMesh = _boundaryMesh->dmMesh();assert(subMesh);
  PetscInt cStart, cEnd;
  PetscErrorCode err = 0;
  err = DMPlexGetHeightStratum(subMesh, 1, &cStart, &cEnd);CHECK_PETSC_ERROR(err);

  const int spaceDim = _quadrature->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();


  topology::Field<topology::SubMesh>& valueField = _parameters->get("value");
  PetscScalar* valueArray = valueField.getLocalArray();

  topology::Field<topology::SubMesh>* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::Field<topology::SubMesh>* rateField = (_dbRate) ? &_parameters->get("rate") : 0;
  topology::Field<topology::SubMesh>* rateTimeField = (_dbRate) ? &_parameters->get("rate time") : 0;
  topology::Field<topology::SubMesh>* changeField = (_dbChange) ? &_parameters->get("change") : 0;
  topology::Field<topology::SubMesh>* changeTimeField = (_dbChange) ? &_parameters->get("change time") : 0;

  PetscScalar* initialArray = (_dbInitial) ? initialField->getLocalArray() : NULL;
  PetscScalar* rateArray = (_dbRate) ? rateField->getLocalArray() : NULL;
  PetscScalar* rateTimeArray = (_dbRate) ? rateTimeField->getLocalArray() : NULL;
  PetscScalar* changeArray = (_dbChange) ? changeField->getLocalArray() : NULL;
  PetscScalar* changeTimeArray = (_dbChange) ? changeTimeField->getLocalArray() : NULL;

  for(PetscInt c = cStart; c < cEnd; ++c) {
    const PetscInt voff = valueField.sectionOffset(c);
    const PetscInt vdof = valueField.sectionDof(c);
    assert(numQuadPts*spaceDim == vdof);
    for (PetscInt d = 0; d < vdof; ++d) {
      valueArray[voff+d] = 0.0;
    } // for

    // Contribution from initial value
    if (_dbInitial) {
      assert(initialField);
      const PetscInt ioff = initialField->sectionOffset(c);
      const PetscInt idof = initialField->sectionDof(c);
      assert(numQuadPts*spaceDim == idof);
      for (PetscInt d = 0; d < idof; ++d) {
        valueArray[voff+d] += initialArray[ioff+d];
      } // for
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateField);
      const PetscInt roff = rateField->sectionOffset(c);
      const PetscInt rdof = rateField->sectionDof(c);
      assert(numQuadPts*spaceDim == rdof);
      assert(rateTimeField);
      const PetscInt rtoff = rateTimeField->sectionOffset(c);
      assert(numQuadPts == rateTimeField->sectionDof(c));

      for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const PylithScalar tRel = t - rateTimeArray[rtoff+iQuad];
        if (tRel > 0.0)  // rate of change integrated over time
          for(int iDim = 0; iDim < spaceDim; ++iDim) {
            valueArray[voff+iQuad*spaceDim+iDim] += rateArray[roff+iQuad*spaceDim+iDim] * tRel;
	  } // for
      } // for
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      assert(changeField);
      const PetscInt coff = changeField->sectionOffset(c);
      const PetscInt cdof = changeField->sectionDof(c);
      assert(numQuadPts*spaceDim == cdof);
      assert(changeTimeField);
      const PetscInt ctoff = changeTimeField->sectionOffset(c);
      assert(numQuadPts == changeTimeField->sectionDof(c));

      for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const PylithScalar tRel = t - changeTimeArray[ctoff+iQuad];
        if (tRel >= 0) { // change in value over time
          PylithScalar scale = 1.0;
          if (_dbTimeHistory) {
            PylithScalar tDim = tRel;
            _getNormalizer().dimensionalize(&tDim, 1, timeScale);
            const int err = _dbTimeHistory->query(&scale, tDim);
            if (err) {
              std::ostringstream msg;
              msg << "Error querying for time '" << tDim 
                  << "' in time history database "
                  << _dbTimeHistory->label() << ".";
              throw std::runtime_error(msg.str());
            } // if
          } // if
          for (int iDim = 0; iDim < spaceDim; ++iDim) {
            valueArray[voff+iQuad*spaceDim+iDim] += changeArray[coff+iQuad*spaceDim+iDim] * scale;
	  } // for
        } // if
      } // for
    } // if
  } // for
  valueField.restoreLocalArray(&valueArray);
  if (_dbInitial) {
    assert(initialField);
    initialField->restoreLocalArray(&initialArray);
  } // if
  if (_dbRate) {
    assert(rateField);
    rateField->restoreLocalArray(&rateArray);
    assert(rateTimeField);
    rateTimeField->restoreLocalArray(&rateTimeArray);
  } // if
  if (_dbChange) {
    assert(changeField);
    changeField->restoreLocalArray(&changeArray);
    assert(changeTimeField);
    changeTimeField->restoreLocalArray(&changeTimeArray);
  } // if
}  // _calculateValue


// End of file 
