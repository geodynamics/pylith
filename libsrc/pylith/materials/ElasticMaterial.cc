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
typedef pylith::topology::Mesh::RealSection RealSection;

typedef pylith::topology::Field<pylith::topology::Mesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::Mesh>::UpdateAddVisitor UpdateAddVisitor;

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
pylith::materials::ElasticMaterial::initialize(
			    const topology::Mesh& mesh,
			    feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  Material::initialize(mesh, quadrature);

  assert(0 != quadrature);
  _numQuadPts = quadrature->numQuadPts();

  if (0 != _dbInitialStress || 0 != _dbInitialStrain) {
    delete _initialFields; 
    _initialFields =
      new topology::Fields<topology::Field<topology::Mesh> >(mesh);
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
  assert(0 != _properties);
  assert(0 != _stateVars);

  PetscSection   propertiesSection = _properties->petscSection();
  Vec            propertiesVec     = _properties->localVector();
  PetscScalar   *propertiesArray;
  PetscInt       dof, off;
  PetscErrorCode err;

  assert(propertiesSection);assert(propertiesVec);
  err = PetscSectionGetDof(propertiesSection, cell, &dof);CHECK_PETSC_ERROR(err);
  err = PetscSectionGetOffset(propertiesSection, cell, &off);CHECK_PETSC_ERROR(err);
  assert(_propertiesCell.size() == dof);
  err = VecGetArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);
  for(PetscInt d = 0; d < dof; ++d) {_propertiesCell[d] = propertiesArray[off+d];}
  err = VecRestoreArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);

  if (hasStateVars()) {
    PetscSection   stateVarsSection = _stateVars->petscSection();
    Vec            stateVarsVec     = _stateVars->localVector();
    PetscScalar   *stateVarsArray;
    PetscInt       dof, off;
    PetscErrorCode err;

    assert(stateVarsSection);assert(stateVarsVec);
    err = PetscSectionGetDof(stateVarsSection, cell, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(stateVarsSection, cell, &off);CHECK_PETSC_ERROR(err);
    assert(_stateVarsCell.size() == dof);
    err = VecGetArray(stateVarsVec, &stateVarsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d) {_stateVarsCell[d] = stateVarsArray[off+d];}
    err = VecRestoreArray(stateVarsVec, &stateVarsArray);CHECK_PETSC_ERROR(err);
  } // if

  _initialStressCell = 0.0;
  _initialStrainCell = 0.0;
  if (0 != _initialFields) {
    if (_initialFields->hasField("initial stress")) {
      PetscSection   stressSection = _initialFields->get("initial stress").petscSection();
      Vec            stressVec     = _initialFields->get("initial stress").localVector();
      PetscScalar   *stressArray;
      PetscInt       dof, off;
      PetscErrorCode err;

      assert(stressSection);assert(stressVec);
      err = PetscSectionGetDof(stressSection, cell, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(stressSection, cell, &off);CHECK_PETSC_ERROR(err);
      assert(_initialStressCell.size() == dof);
      err = VecGetArray(stressVec, &stressArray);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {_initialStressCell[d] = stressArray[off+d];}
      err = VecRestoreArray(stressVec, &stressArray);CHECK_PETSC_ERROR(err);
    } // if
    if (_initialFields->hasField("initial strain")) {
      PetscSection   strainSection = _initialFields->get("initial strain").petscSection();
      Vec            strainVec     = _initialFields->get("initial strain").localVector();
      PetscScalar   *strainArray;
      PetscInt       dof, off;
      PetscErrorCode err;

      assert(strainSection);assert(strainVec);
      err = PetscSectionGetDof(strainSection, cell, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(strainSection, cell, &off);CHECK_PETSC_ERROR(err);
      assert(_initialStrainCell.size() == dof);
      err = VecGetArray(strainVec, &strainArray);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d) {_initialStrainCell[d] = strainArray[off+d];}
      err = VecRestoreArray(strainVec, &strainArray);CHECK_PETSC_ERROR(err);
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
pylith::materials::ElasticMaterial::calcDerivElastic(
					       const scalar_array& totalStrain)
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
pylith::materials::ElasticMaterial::updateStateVars(
					      const scalar_array& totalStrain,
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
  
  const ALE::Obj<RealSection>& stateVarsSection = _stateVars->section();
  assert(!stateVarsSection.isNull());  
  stateVarsSection->updatePoint(cell, &_stateVarsCell[0]);
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
  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Setup field if necessary.
  PetscSection fieldSection = PETSC_NULL;
  PetscVec fieldVec = PETSC_NULL;
  PetscScalar* fieldArray = PETSC_NULL;
  PetscInt dof, off;
  if (field) {
    const int fiberDim = 1*numQuadPts;
    fieldSection = field->petscSection();
    bool useCurrentField = false;
    if (fieldSection) {
      // check fiber dimension
      int fiberDimCurrentLocal = 0;
      int fiberDimCurrent = 0;
      if (numCells > 0) {
	err = PetscSectionGetDof(fieldSection, cells[0], &fiberDimCurrentLocal);CHECK_PETSC_ERROR(err);
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
      fieldSection = field->petscSection();
      logger.stagePop();
    } // if
    assert(fieldSection);
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    err = VecGetArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
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
      err = PetscSectionGetDof(fieldSection, cell, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(fieldSection, cell, &off);CHECK_PETSC_ERROR(err);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
	fieldArray[off+iQuad] = dtStableCell[iQuad];
      } // for
    } // if
  } // for
  if (field) {
    err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
  } // if

  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  
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
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Setup field if necessary.
  ALE::Obj<RealSection> fieldSection;
  if (field) {
    const int fiberDim = 1*numQuadPts;
    fieldSection = field->section();
    bool useCurrentField = false;
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int fiberDimCurrentLocal = (cells->size() > 0) ? fieldSection->getFiberDimension(*cells->begin()) : 0;
      int fiberDimCurrent = 0;
      MPI_Allreduce((void *) &fiberDimCurrentLocal, 
		    (void *) &fiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      field->newSection(cells, fiberDim);
      field->allocate();
      fieldSection = field->section();
      logger.stagePop();
    } // if
    assert(!fieldSection.isNull());
    field->label("stable_dt_explicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
  } // if

  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();
  scalar_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  scalar_array dtStableCell(numQuadPts);
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    retrievePropsAndVars(*c_iter);

    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
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
      assert(!fieldSection.isNull());
      assert(numQuadPts == fieldSection->getFiberDimension(*c_iter));
      fieldSection->updatePoint(*c_iter, &dtStableCell[0]);
    } // if
  } // for
  
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

    const ALE::Obj<RealSection>& fieldSection = field->section();
    // Get cells associated with material
    const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
    assert(!sieveMesh.isNull());
    const ALE::Obj<SieveMesh::label_sequence>& cells = sieveMesh->getLabelStratum("material-id", id());
    assert(!cells.isNull());
    const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
    const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
    
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (!fieldSection.isNull()) {
      // check fiber dimension
      const int fiberDimCurrentLocal = (cells->size() > 0) ? fieldSection->getFiberDimension(*cells->begin()) : 0;
      int fiberDimCurrent = 0;
      MPI_Allreduce((void *) &fiberDimCurrentLocal, 
		    (void *) &fiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      field->newSection(cells, fiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    assert(!fieldSection.isNull());
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);

    scalar_array dtStableCell(numQuadPts);
    dtStableCell = PYLITH_MAXSCALAR;
    for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
	 ++c_iter) {
      assert(numQuadPts == fieldSection->getFiberDimension(*c_iter));
      fieldSection->updatePoint(*c_iter, &dtStableCell[0]);
    } // for
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
pylith::materials::ElasticMaterial::_initializeInitialStress(
			 const topology::Mesh& mesh,
			 feassemble::Quadrature<topology::Mesh>* quadrature)
{ // _initializeInitialStress
  if (0 == _dbInitialStress)
    return;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("MaterialsFields");

  assert(0 != _initialFields);
  _initialFields->add("initial stress", "initial_stress");
  topology::Field<topology::Mesh>& initialStress = 
    _initialFields->get("initial stress");

  assert(0 != _dbInitialStress);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  // Get cells associated with material
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
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
  Vec          initialStressVec     = initialStress.localVector();

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
  
  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
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

    err = DMComplexVecSetClosure(dmMesh, initialStressSection, initialStressVec, cell, &stressCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
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

  assert(0 != _initialFields);
  _initialFields->add("initial strain", "initial_strain");
  topology::Field<topology::Mesh>& initialStrain = 
    _initialFields->get("initial strain");

  assert(0 != _dbInitialStrain);
  assert(0 != quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  // Get cells associated with material
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
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
  Vec          initialStrainVec     = initialStrain.localVector();

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
  
  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
    
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
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

    err = DMComplexVecSetClosure(dmMesh, initialStrainSection, initialStrainVec, cell, &strainCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
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
