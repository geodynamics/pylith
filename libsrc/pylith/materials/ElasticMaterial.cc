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

#include "ElasticMaterial.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial::ElasticMaterial(const int dimension,
						    const int tensorSize,
						    const int numElasticConsts,
						    const Metadata& metadata) :
  Material(dimension, tensorSize, metadata),
  _dbInitialStress(0),
  _dbInitialStrain(0),
  _initialFields(0),
  _numQuadPts(0),
  _numElasticConsts(numElasticConsts)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial::~ElasticMaterial(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::ElasticMaterial::deallocate(void)
{ // deallocate
  Material::deallocate();
  delete _initialFields; _initialFields = 0;

  _dbInitialStress = 0; // :TODO: Use shared pointer.
  _dbInitialStrain = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize material by getting physical property parameters from
// database.
void
pylith::materials::ElasticMaterial::initialize(const topology::Mesh& mesh,
					       feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  Material::initialize(mesh, quadrature);

  assert(0 != quadrature);
  _numQuadPts = quadrature->numQuadPts();

  if (_dbInitialStress || _dbInitialStrain) {
    delete _initialFields; 
    _initialFields = new topology::Fields<topology::Field<topology::Mesh> >(mesh);
  } // if
  _initializeInitialStress(mesh, quadrature);
  _initializeInitialStrain(mesh, quadrature);
  _allocateCellArrays();
} // initialize

// ----------------------------------------------------------------------
// Retrieve parameters for physical properties and state variables for cell.
void
pylith::materials::ElasticMaterial::retrievePropsAndVars(const int cell)
{ // retrievePropsAndVars
  assert(_properties);
  assert(_stateVars);

  const int propertiesSize = _numQuadPts*_numPropsQuadPt;
  const int stateVarsSize = _numQuadPts*_numVarsQuadPt;
  assert(_propertiesCell.size() == propertiesSize);
  assert(_stateVarsCell.size() == stateVarsSize);

  topology::VecVisitorMesh propertiesVisitor(*_properties);
  PetscScalar* propertiesArray = propertiesVisitor.localArray();
  const PetscInt poff = propertiesVisitor.sectionOffset(cell);
  assert(propertiesSize == propertiesVisitor.sectionDof(cell));
  for(PetscInt d = 0; d < propertiesSize; ++d) {
    _propertiesCell[d] = propertiesArray[poff+d];
  } // for

  if (hasStateVars()) {
    topology::VecVisitorMesh stateVarsVisitor(*_stateVars);
    PetscScalar* stateVarsArray = stateVarsVisitor.localArray();
    const PetscInt soff = stateVarsVisitor.sectionOffset(cell);
    assert(stateVarsSize == stateVarsVisitor.sectionDof(cell));
    for(PetscInt d = 0; d < stateVarsSize; ++d) {
      _stateVarsCell[d] = stateVarsArray[soff+d];
    } // for
  } // if

  _initialStressCell = 0.0;
  _initialStrainCell = 0.0;
  if (_initialFields) {
    if (_initialFields->hasField("initial stress")) {
      const int stressSize = _numQuadPts*_tensorSize;
      assert(_initialStressCell.size() == stressSize);
      topology::Field<topology::Mesh>& stressField = _initialFields->get("initial stress");
      topology::VecVisitorMesh stressVisitor(stressField);
      PetscScalar* stressArray = stressVisitor.localArray();
      const PetscInt ioff = stressVisitor.sectionOffset(cell);
      assert(stressSize == stressVisitor.sectionDof(cell));
      for(PetscInt d = 0; d < stressSize; ++d) {
	_initialStressCell[d] = stressArray[ioff+d];
      } // for
    } // if
    if (_initialFields->hasField("initial strain")) {
      const int strainSize = _numQuadPts*_tensorSize;
      assert(_initialStrainCell.size() == strainSize);
      topology::Field<topology::Mesh>& strainField = _initialFields->get("initial strain");
      topology::VecVisitorMesh strainVisitor(strainField);
      PetscScalar* strainArray = strainVisitor.localArray();
      const PetscInt ioff = strainVisitor.sectionOffset(cell);
      assert(strainSize == strainVisitor.sectionDof(cell));
      for(PetscInt d = 0; d < strainSize; ++d) {
	_initialStrainCell[d] = strainArray[ioff+d];
      } // for
    } // if
  } // if
} // retrievePropsAndVars

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcStress(const scalar_array& totalStrain,
					       const bool computeStateVars)
{ // calcStress
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_stressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stressCell[iQuad*_tensorSize], _tensorSize,
		&_propertiesCell[iQuad*numPropsQuadPt], numPropsQuadPt,
		&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		&totalStrain[iQuad*_tensorSize], _tensorSize, 
		&_initialStressCell[iQuad*_tensorSize], _tensorSize,
		&_initialStrainCell[iQuad*_tensorSize], _tensorSize,
		computeStateVars);

  return _stressCell;
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcDerivElastic(const scalar_array& totalStrain)
{ // calcDerivElastic
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_elasticConstsCell.size() == numQuadPts*_numElasticConsts);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConstsCell[iQuad*_numElasticConsts], 
		       _numElasticConsts,
		       &_propertiesCell[iQuad*numPropsQuadPt], 
		       numPropsQuadPt, 
		       &_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		       &totalStrain[iQuad*_tensorSize], _tensorSize,
		       &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		       &_initialStrainCell[iQuad*_tensorSize], _tensorSize);

  return _elasticConstsCell;
} // calcDerivElastic

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::materials::ElasticMaterial::updateStateVars(const scalar_array& totalStrain,
						    const int cell)
{ // updateStateVars
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);
  assert(totalStrain.size() == numQuadPts*_tensorSize);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _updateStateVars(&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		     &_propertiesCell[iQuad*numPropsQuadPt], 
		     numPropsQuadPt,
		     &totalStrain[iQuad*_tensorSize], _tensorSize,
		     &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		     &_initialStrainCell[iQuad*_tensorSize], _tensorSize);
  
  topology::VecVisitorMesh stateVarsVisitor(*_stateVars);
  PetscScalar* stateVarsArray = stateVarsVisitor.localArray();
  const PetscInt soff = stateVarsVisitor.sectionOffset(cell);
  const int stateVarsSize = numQuadPts*numVarsQuadPt;
  assert(stateVarsSize == stateVarsVisitor.sectionDof(cell));
  for (PetscInt d = 0; d < stateVarsSize; ++d) {
    stateVarsArray[soff+d] = _stateVarsCell[d];
  } // for
} // updateStateVars

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticMaterial::stableTimeStepImplicit(const topology::Mesh& mesh,
							   topology::Field<topology::Mesh>* field)
{ // stableTimeStepImplicit
  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_elasticConstsCell.size() == numQuadPts*_numElasticConsts);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", id());
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();
 
  // Setup field if necessary.
  topology::VecVisitorMesh* fieldVisitor = (field) ? new topology::VecVisitorMesh(*field) : 0;
  PetscScalar* fieldArray = NULL;
  if (field) {
    assert(fieldVisitor);
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (fieldVisitor->petscSection()) {
      // check fiber dimension
      int fiberDimCurrentLocal = 0;
      int fiberDimCurrent = 0;
      if (numCells > 0) {
	fiberDimCurrentLocal = fieldVisitor->sectionDof(cells[0]);
      } // if
      MPI_Allreduce((void *) &fiberDimCurrentLocal, 
		    (void *) &fiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      int_array cellsTmp(cells, numCells);
      field->newSection(cellsTmp, fiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    fieldArray = fieldVisitor->localArray();
  } // if

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  scalar_array dtStableCell(numQuadPts);
  for (PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    retrievePropsAndVars(cell);
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar dt = 
        _stableTimeStepImplicit(&_propertiesCell[iQuad*numPropsQuadPt],
                                numPropsQuadPt,
                                &_stateVarsCell[iQuad*numVarsQuadPt],
                                numVarsQuadPt);
      dtStableCell[iQuad] = dt;
      if (dt < dtStable) {
        dtStable = dt;
      } // if
    } // for
    if (field) {
      assert(fieldVisitor);
      const PetscInt off = fieldVisitor->sectionOffset(cell);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	fieldArray[off+iQuad] = dtStableCell[iQuad];
      } // for
    } // if
  } // for
  delete fieldVisitor; fieldVisitor = 0;

  return dtStable;
} // stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::ElasticMaterial::stableTimeStepExplicit(const topology::Mesh& mesh,
							   feassemble::Quadrature<topology::Mesh>* quadrature,
							   topology::Field<topology::Mesh>* field)
{ // stableTimeStepImplicit
  assert(quadrature);

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int tensorSize = _tensorSize;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == numQuadPts*numPropsQuadPt);
  assert(_stateVarsCell.size() == numQuadPts*numVarsQuadPt);
  assert(_elasticConstsCell.size() == numQuadPts*_numElasticConsts);
  assert(_initialStressCell.size() == numQuadPts*_tensorSize);
  assert(_initialStrainCell.size() == numQuadPts*_tensorSize);

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::StratumIS materialIS(dmMesh, "material-id", id());
  const PetscInt* cells = materialIS.points();
  const PetscInt numCells = materialIS.size();

  // Setup field if necessary.
  topology::VecVisitorMesh* fieldVisitor = (field) ? new topology::VecVisitorMesh(*field) : 0;
  PetscScalar *fieldArray = NULL;
  if (field) {
    assert(fieldVisitor);
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (fieldVisitor->petscSection()) {
      // check fiber dimension
      PetscInt fiberDimCurrentLocal = 0;
      if (numCells > 0) {
	fiberDimCurrentLocal = fieldVisitor->sectionDof(cells[0]);
      } // if
      PetscInt fiberDimCurrent = 0;
      MPI_Allreduce(&fiberDimCurrentLocal, &fiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      field->newSection(cells, numCells, fiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    field->label("stable_dt_explicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    fieldArray = fieldVisitor->localArray();
  } // if

  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  scalar_array coordinatesCell(numBasis*spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar *coordsArray = NULL;
  PetscInt coordsSize = 0;

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  scalar_array dtStableCell(numQuadPts);
  for (PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    retrievePropsAndVars(cell);

    coordsVisitor.getClosure(&coordsArray, &coordsSize, cell);
    for(PetscInt i = 0; i < coordsSize; ++i) {
      coordinatesCell[i] = coordsArray[i];
    } // for
    coordsVisitor.restoreClosure(&coordsArray, &coordsSize, cell);
    const double minCellWidth = quadrature->minCellWidth(coordinatesCell);

    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const PylithScalar dt = 
	_stableTimeStepExplicit(&_propertiesCell[iQuad*numPropsQuadPt],
				numPropsQuadPt,
				&_stateVarsCell[iQuad*numVarsQuadPt],
				numVarsQuadPt,
				minCellWidth);
      dtStableCell[iQuad] = dt;
      if (dt < dtStable) {
        dtStable = dt;
      } // if
    } // for
    if (field) {
      assert(fieldVisitor);
      const PetscInt off = fieldVisitor->sectionOffset(cell);
      assert(numQuadPts == fieldVisitor->sectionDof(cell));
      for (PetscInt d = 0; d < numQuadPts; ++d) {
        fieldArray[off+d] = dtStableCell[d];
      } // for
    } // if
  } // for
  delete fieldVisitor; fieldVisitor = 0;

  return dtStable;
} // stableTimeStepExplicit

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration (return large value).
PylithScalar
pylith::materials::ElasticMaterial::_stableTimeStepImplicitMax(const topology::Mesh& mesh,
							       topology::Field<topology::Mesh>* field)
{ // _stableTimeStepImplicitMax
  const PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  
  if (field) {
    const int numQuadPts = _numQuadPts;

    PetscSection fieldSection = field->petscSection();
    PetscVec fieldVec = field->localVector();
    PetscScalar *fieldArray;
    // Get cells associated with material
    PetscDM dmMesh = mesh.dmMesh();
    assert(dmMesh);
    PetscIS cellIS;
    const PetscInt *cells;
    PetscInt numCells;
    PetscErrorCode err = 0;

    err = DMPlexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
    err = ISGetLocalSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (fieldSection) {
      // check fiber dimension
      PetscInt fiberDimCurrentLocal = 0;
      if (numCells > 0) {err = PetscSectionGetDof(fieldSection, cells[0], &fiberDimCurrentLocal);CHECK_PETSC_ERROR(err);}
      PetscInt fiberDimCurrent = 0;
      MPI_Allreduce(&fiberDimCurrentLocal, &fiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      field->newSection(cells, numCells, fiberDim);
      field->allocate();
      fieldSection = field->petscSection();
      fieldVec     = field->localVector();
      logger.stagePop();
    } // if
    assert(fieldSection);
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);

    scalar_array dtStableCell(numQuadPts);
    dtStableCell = PYLITH_MAXSCALAR;
    err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
    for (PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];
      PetscInt dof, off;

      assert(fieldSection);
      err = PetscSectionGetDof(fieldSection, cell, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(fieldSection, cell, &off);CHECK_PETSC_ERROR(err);
      assert(numQuadPts == dof);
      for (PetscInt d = 0; d < dof; ++d) {
        fieldArray[off+d] = dtStableCell[d];
      }
    } // for
    err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
    err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  } // if
  
  return dtStable;
} // _stableTimeStepImplicitMax

// ----------------------------------------------------------------------
// Allocate cell arrays.
void
pylith::materials::ElasticMaterial::_allocateCellArrays(void)
{ // _allocateCellArrays
  const int numQuadPts = _numQuadPts;
  const int tensorSize = _tensorSize;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  const int numElasticConsts = _numElasticConsts;

  _propertiesCell.resize(numQuadPts * numPropsQuadPt);
  _stateVarsCell.resize(numQuadPts * numVarsQuadPt);
  _initialStressCell.resize(numQuadPts * tensorSize);
  _initialStrainCell.resize(numQuadPts * tensorSize);
  _densityCell.resize(numQuadPts);
  _stressCell.resize(numQuadPts * tensorSize);
  _elasticConstsCell.resize(numQuadPts * numElasticConsts);
} // _allocateCellArrays

// ----------------------------------------------------------------------
// Initialize initial stress field.
void
pylith::materials::ElasticMaterial::_initializeInitialStress(const topology::Mesh& mesh,
							     feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStress
  if (0 == _dbInitialStress)
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("MaterialsFields");

  assert(_initialFields);
  _initialFields->add("initial stress", "initial_stress");
  topology::Field<topology::Mesh>& initialStress = 
    _initialFields->get("initial stress");

  assert(_dbInitialStress);
  assert(quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  // Get cells associated with material
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  PetscDM dmMesh = mesh.dmMesh();
  PetscIS cellIS;
  const PetscInt *cells;
  PetscInt numCells;
  PetscErrorCode err = 0;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array stressCell(numQuadPts*tensorSize);

  // Create field to hold initial stress state.
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  int_array cellsTmp(cells, numCells);
  initialStress.newSection(cellsTmp, fiberDim);
  initialStress.allocate();
  initialStress.zero();
  PetscSection initialStressSection = initialStress.petscSection();
  PetscVec initialStressVec = initialStress.localVector();

  // Setup databases for querying
  _dbInitialStress->open();
  switch (dimension())
    { // switch
    case 1: {
      const char* stressDBValues[] = { "stress" };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 1
    case 2 : {
      const char* stressDBValues[] = { 
	"stress-xx",
	"stress-yy",
	"stress-xy",
      };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 2
    case 3 : {
      const char* stressDBValues[] = { 
	"stress-xx",
	"stress-yy",
	"stress-zz",
	"stress-xy",
	"stress-yz",
	"stress-xz",
      };
      _dbInitialStress->queryVals(stressDBValues, tensorSize);
      break;
    } // case 3
    default :
      std::cerr << "Bad dimension '" << dimension() << "'." << std::endl;
      assert(0);
      throw std::logic_error("Unknown dimension in elastic material.");
    } // switch
  
  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    PetscScalar *coords;
    PetscInt     coordsSize;
    err = DMPlexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    quadrature->computeGeometry(coordinatesCell, cell);
    err = DMPlexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, iCoord=0, iStress=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, iCoord+=spaceDim, iStress+=tensorSize) {
      int err = _dbInitialStress->query(&stressCell[iStress], tensorSize,
					&quadPtsGlobal[iCoord], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial stress at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iCoord+i];
	msg << ") in material " << label() << "\n"
	    << "using spatial database '" << _dbInitialStress->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    // Nondimensionalize stress
    _normalizer->nondimensionalize(&stressCell[0], stressCell.size(), 
				   pressureScale);

    err = DMPlexVecSetClosure(dmMesh, initialStressSection, initialStressVec, cell, &stressCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
  } // for

  // Close databases
  _dbInitialStress->close();

  logger.stagePop();
} // _initializeInitialStress

// ----------------------------------------------------------------------
// Initialize initial strain field.
void
pylith::materials::ElasticMaterial::_initializeInitialStrain(
			 const topology::Mesh& mesh,
			 feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStrain
  if (0 == _dbInitialStrain)
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("MaterialsFields");

  assert(_initialFields);
  _initialFields->add("initial strain", "initial_strain");
  topology::Field<topology::Mesh>& initialStrain = 
    _initialFields->get("initial strain");

  assert(_dbInitialStrain);
  assert(quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  // Get cells associated with material
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  PetscDM dmMesh = mesh.dmMesh();
  PetscIS cellIS;
  const PetscInt *cells;
  PetscInt numCells;
  PetscErrorCode err = 0;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  // Create arrays for querying
  const int tensorSize = _tensorSize;
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);

  // Create field to hold initial strain state.
  const int fiberDim = numQuadPts * tensorSize;
  assert(fiberDim > 0);
  int_array cellsTmp(cells, numCells);
  initialStrain.newSection(cellsTmp, fiberDim);
  initialStrain.allocate();
  initialStrain.zero();
  PetscSection initialStrainSection = initialStrain.petscSection();
  PetscVec initialStrainVec = initialStrain.localVector();

  // Setup databases for querying
  _dbInitialStrain->open();
  switch (dimension())
    { // switch
    case 1: {
      const char* strainDBValues[] = { "strain" };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 1
    case 2 : {
      const char* strainDBValues[] = { 
	"strain-xx",
	"strain-yy",
	"strain-xy",
      };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 2
    case 3 : {
      const char* strainDBValues[] = { 
	"strain-xx",
	"strain-yy",
	"strain-zz",
	"strain-xy",
	"strain-yz",
	"strain-xz",
      };
      _dbInitialStrain->queryVals(strainDBValues, tensorSize);
      break;
    } // case 3
    default :
      std::cerr << "Bad dimension '" << dimension() << "'." << std::endl;
      assert(0);
      throw std::logic_error("Unknown dimension in elastic material.");
    } // switch
  
  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
    
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    PetscScalar *coords;
    PetscInt     coordsSize;
    err = DMPlexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    quadrature->computeGeometry(coordinatesCell, cell);
    err = DMPlexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);
    
    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, iCoord=0, iStrain=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, iCoord+=spaceDim, iStrain+=tensorSize) {
      int err = _dbInitialStrain->query(&strainCell[iStrain], tensorSize,
					&quadPtsGlobal[iCoord], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial strain at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iCoord+i];
	msg << ") in material " << label() << "\n"
	    << "using spatial database '" << _dbInitialStrain->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    err = DMPlexVecSetClosure(dmMesh, initialStrainSection, initialStrainVec, cell, &strainCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
  } // for

  // Close databases
  _dbInitialStrain->close();

  logger.stagePop();
} // _initializeInitialStrain

// ----------------------------------------------------------------------
// Update stateVars (for next time step).
void
pylith::materials::ElasticMaterial::_updateStateVars(
					    PylithScalar* const stateVars,
					    const int numStateVars,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* totalStrain,
					    const int strainSize,
					    const PylithScalar* initialStress,
					    const int initialStressSize,
					    const PylithScalar* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVars
} // _updateStateVars


// End of file 
