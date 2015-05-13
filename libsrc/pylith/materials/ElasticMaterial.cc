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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ElasticMaterial.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

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
  PYLITH_METHOD_BEGIN;

  Material::deallocate();
  delete _initialFields; _initialFields = 0;

  _dbInitialStress = 0; // :TODO: Use shared pointer.
  _dbInitialStrain = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize material by getting physical property parameters from
// database.
void
pylith::materials::ElasticMaterial::initialize(const topology::Mesh& mesh,
					       feassemble::Quadrature* quadrature)
{ // initialize
  PYLITH_METHOD_BEGIN;

  Material::initialize(mesh, quadrature);

  assert(0 != quadrature);
  _numQuadPts = quadrature->numQuadPts();

  if (_dbInitialStress || _dbInitialStrain) {
    delete _initialFields; 
    _initialFields = new topology::Fields(mesh);assert(_initialFields);
  } // if
  _initializeInitialStress(mesh, quadrature);
  _initializeInitialStrain(mesh, quadrature);
  _allocateCellArrays();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Retrieve parameters for physical properties and state variables for cell.
void
pylith::materials::ElasticMaterial::retrievePropsAndVars(const int cell)
{ // retrievePropsAndVars
  PYLITH_METHOD_BEGIN;

  assert(_properties);
  assert(_stateVars);

  const int propertiesSize = _numQuadPts*_numPropsQuadPt;
  const int stateVarsSize = _numQuadPts*_numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(propertiesSize));
  assert(_stateVarsCell.size() == size_t(stateVarsSize));

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
      assert(_initialStressCell.size() == size_t(stressSize));
      topology::Field& stressField = _initialFields->get("initial stress");
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
      assert(_initialStrainCell.size() == size_t(strainSize));
      topology::Field& strainField = _initialFields->get("initial strain");
      topology::VecVisitorMesh strainVisitor(strainField);
      PetscScalar* strainArray = strainVisitor.localArray();
      const PetscInt ioff = strainVisitor.sectionOffset(cell);
      assert(strainSize == strainVisitor.sectionDof(cell));
      for(PetscInt d = 0; d < strainSize; ++d) {
	_initialStrainCell[d] = strainArray[ioff+d];
      } // for
    } // if
  } // if

  PYLITH_METHOD_END;
} // retrievePropsAndVars

// ----------------------------------------------------------------------
// Compute stress tensor for cell at quadrature points.
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcStress(const scalar_array& totalStrain,
					       const bool computeStateVars)
{ // calcStress
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(numQuadPts*numPropsQuadPt));
  assert(_stateVarsCell.size() == size_t(numQuadPts*numVarsQuadPt));
  assert(_stressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStrainCell.size() == size_t(numQuadPts*_tensorSize));
  assert(totalStrain.size() == size_t(numQuadPts*_tensorSize));

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcStress(&_stressCell[iQuad*_tensorSize], _tensorSize,
		&_propertiesCell[iQuad*numPropsQuadPt], numPropsQuadPt,
		&_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		&totalStrain[iQuad*_tensorSize], _tensorSize, 
		&_initialStressCell[iQuad*_tensorSize], _tensorSize,
		&_initialStrainCell[iQuad*_tensorSize], _tensorSize,
		computeStateVars);

  PYLITH_METHOD_RETURN(_stressCell);
} // calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix for cell at quadrature points.
const pylith::scalar_array&
pylith::materials::ElasticMaterial::calcDerivElastic(const scalar_array& totalStrain)
{ // calcDerivElastic
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(numQuadPts*numPropsQuadPt));
  assert(_stateVarsCell.size() == size_t(numQuadPts*numVarsQuadPt));
  assert(_elasticConstsCell.size() == size_t(numQuadPts*_numElasticConsts));
  assert(_initialStressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStrainCell.size() == size_t(numQuadPts*_tensorSize));
  assert(totalStrain.size() == size_t(numQuadPts*_tensorSize));

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    _calcElasticConsts(&_elasticConstsCell[iQuad*_numElasticConsts], 
		       _numElasticConsts,
		       &_propertiesCell[iQuad*numPropsQuadPt], 
		       numPropsQuadPt, 
		       &_stateVarsCell[iQuad*numVarsQuadPt], numVarsQuadPt,
		       &totalStrain[iQuad*_tensorSize], _tensorSize,
		       &_initialStressCell[iQuad*_tensorSize], _tensorSize,
		       &_initialStrainCell[iQuad*_tensorSize], _tensorSize);

  PYLITH_METHOD_RETURN(_elasticConstsCell);
} // calcDerivElastic

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::materials::ElasticMaterial::updateStateVars(const scalar_array& totalStrain,
						    const int cell)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(numQuadPts*numPropsQuadPt));
  assert(_stateVarsCell.size() == size_t(numQuadPts*numVarsQuadPt));
  assert(_initialStressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStrainCell.size() == size_t(numQuadPts*_tensorSize));
  assert(totalStrain.size() == size_t(numQuadPts*_tensorSize));

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

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticMaterial::stableTimeStepImplicit(const topology::Mesh& mesh,
							   topology::Field* field)
{ // stableTimeStepImplicit
  PYLITH_METHOD_BEGIN;

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(numQuadPts*numPropsQuadPt));
  assert(_stateVarsCell.size() == size_t(numQuadPts*numVarsQuadPt));
  assert(_elasticConstsCell.size() == size_t(numQuadPts*_numElasticConsts));
  assert(_initialStressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStrainCell.size() == size_t(numQuadPts*_tensorSize));

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();
 
  // Setup field if necessary.
  topology::VecVisitorMesh* fieldVisitor = NULL;
  PetscScalar* fieldArray = NULL;
  if (field) {
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (field->hasSection()) {
      // check fiber dimension
      int fiberDimCurrentLocal = 0;
      int fiberDimCurrent = 0;
      if (numCells > 0) {
	topology::VecVisitorMesh fieldVisitor(*field);
	fiberDimCurrentLocal = fieldVisitor.sectionDof(cells[0]);
      } // if
      MPI_Allreduce((void *) &fiberDimCurrentLocal, 
		    (void *) &fiberDimCurrent, 1, 
		    MPI_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      int_array cellsTmp(cells, numCells);
      field->newSection(cellsTmp, fiberDim);
      field->allocate();
    } // if
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    fieldVisitor = new topology::VecVisitorMesh(*field);assert(fieldVisitor);
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

  assert(dtStable > 0.0);

  PYLITH_METHOD_RETURN(dtStable);
} // stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::ElasticMaterial::stableTimeStepExplicit(const topology::Mesh& mesh,
							   feassemble::Quadrature* quadrature,
							   topology::Field* field)
{ // stableTimeStepImplicit
  PYLITH_METHOD_BEGIN;

  assert(quadrature);

  const int numQuadPts = _numQuadPts;
  const int numPropsQuadPt = _numPropsQuadPt;
  const int numVarsQuadPt = _numVarsQuadPt;
  assert(_propertiesCell.size() == size_t(numQuadPts*numPropsQuadPt));
  assert(_stateVarsCell.size() == size_t(numQuadPts*numVarsQuadPt));
  assert(_elasticConstsCell.size() == size_t(numQuadPts*_numElasticConsts));
  assert(_initialStressCell.size() == size_t(numQuadPts*_tensorSize));
  assert(_initialStrainCell.size() == size_t(numQuadPts*_tensorSize));

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup field if necessary.
  topology::VecVisitorMesh* fieldVisitor = NULL;
  PetscScalar *fieldArray = NULL;
  if (field) {
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (field->hasSection()) {
      // check fiber dimension
      PetscInt fiberDimCurrentLocal = 0;
      if (numCells > 0) {
	topology::VecVisitorMesh fieldVisitor(*field);
	fiberDimCurrentLocal = fieldVisitor.sectionDof(cells[0]);
      } // if
      PetscInt fiberDimCurrent = 0;
      MPI_Allreduce(&fiberDimCurrentLocal, &fiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, numCells, fiberDim);
      field->allocate();
    } // if
    field->label("stable_dt_explicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    fieldVisitor = new topology::VecVisitorMesh(*field);assert(fieldVisitor);
    fieldArray = fieldVisitor->localArray();
  } // if

  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  scalar_array dtStableCell(numQuadPts);
  for (PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    retrievePropsAndVars(cell);

    coordsVisitor.getClosure(&coordsCell, cell);
    const PylithScalar minCellWidth = quadrature->minCellWidth(&coordsCell[0], numBasis, spaceDim);
  assert(minCellWidth > 0.0);

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

  assert(dtStable > 0.0);

  PYLITH_METHOD_RETURN(dtStable);
} // stableTimeStepExplicit

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration (return large value).
PylithScalar
pylith::materials::ElasticMaterial::_stableTimeStepImplicitMax(const topology::Mesh& mesh,
							       topology::Field* field)
{ // _stableTimeStepImplicitMax
  PYLITH_METHOD_BEGIN;

  const PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;
  
  if (field) {
    const int numQuadPts = _numQuadPts;

    // Get cells associated with material
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);

    assert(_materialIS);
    const PetscInt* cells = _materialIS->points();
    const PetscInt numCells = _materialIS->size();
    
    // Setup field if necessary.
    const int fiberDim = 1*numQuadPts;
    bool useCurrentField = false;
    if (field->hasSection()) {
      // check fiber dimension
      PetscInt fiberDimCurrentLocal = 0;
      if (numCells > 0) {
	topology::VecVisitorMesh fieldVisitor(*field);
	fiberDimCurrentLocal = fieldVisitor.sectionDof(cells[0]);
      } // if
      PetscInt fiberDimCurrent = 0;
      MPI_Allreduce(&fiberDimCurrentLocal, &fiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(fiberDimCurrent > 0);
      useCurrentField = fiberDim == fiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, numCells, fiberDim);
      field->allocate();
    } // if
    field->label("stable_dt_implicit");
    assert(_normalizer);
    field->scale(_normalizer->timeScale());
    field->vectorFieldType(topology::FieldBase::MULTI_SCALAR);
    topology::VecVisitorMesh fieldVisitor(*field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    scalar_array dtStableCell(numQuadPts);
    dtStableCell = PYLITH_MAXSCALAR;
    for (PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];

      const PetscInt off = fieldVisitor.sectionOffset(cell);
      assert(numQuadPts == fieldVisitor.sectionDof(cell));
      for (PetscInt d = 0; d < numQuadPts; ++d) {
        fieldArray[off+d] = dtStableCell[d];
      } // for
    } // for
  } // if
  
  PYLITH_METHOD_RETURN(dtStable);
} // _stableTimeStepImplicitMax

// ----------------------------------------------------------------------
// Allocate cell arrays.
void
pylith::materials::ElasticMaterial::_allocateCellArrays(void)
{ // _allocateCellArrays
  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_END;
} // _allocateCellArrays

// ----------------------------------------------------------------------
// Initialize initial stress field.
void
pylith::materials::ElasticMaterial::_initializeInitialStress(const topology::Mesh& mesh,
							     feassemble::Quadrature* quadrature)
{ // _initializeInitialStress
  PYLITH_METHOD_BEGIN;

  if (!_dbInitialStress)
    PYLITH_METHOD_END;

  assert(_initialFields);
  _initialFields->add("initial stress", "initial_stress");
  topology::Field& initialStress = _initialFields->get("initial stress");

  assert(_dbInitialStress);
  assert(quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

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
  initialStress.zeroAll();
  topology::VecVisitorMesh stressVisitor(initialStress);

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
      std::ostringstream msg;
      msg << "Bad dimension '" << dimension() << "' in elastic material." << std::endl;
      throw std::logic_error(msg.str());
    } // switch
  
  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);
    
    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, iCoord=0, iStress=0; iQuadPt < numQuadPts; ++iQuadPt, iCoord+=spaceDim, iStress+=tensorSize) {
      int err = _dbInitialStress->query(&stressCell[iStress], tensorSize,
					&quadPtsGlobal[iCoord], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial stress at (";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iCoord+i];
	msg << ") in material '" << label() << "' using spatial database '" << _dbInitialStress->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    // Nondimensionalize stress
    _normalizer->nondimensionalize(&stressCell[0], stressCell.size(), 
				   pressureScale);

    stressVisitor.setClosure(&stressCell[0], stressCell.size(), cell, ADD_VALUES);
  } // for

  // Close databases
  _dbInitialStress->close();

  PYLITH_METHOD_END;
} // _initializeInitialStress

// ----------------------------------------------------------------------
// Initialize initial strain field.
void
pylith::materials::ElasticMaterial::_initializeInitialStrain(const topology::Mesh& mesh,
							     feassemble::Quadrature* quadrature)
{ // _initializeInitialStrain
  PYLITH_METHOD_BEGIN;

  if (!_dbInitialStrain)
    PYLITH_METHOD_END;

  assert(_initialFields);
  _initialFields->add("initial strain", "initial_strain");
  topology::Field& initialStrain = _initialFields->get("initial strain");

  assert(_dbInitialStrain);
  assert(quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int spaceDim = quadrature->spaceDim();
  const int numBasis = quadrature->numBasis();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

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
  initialStrain.zeroAll();
  topology::VecVisitorMesh strainVisitor(initialStrain);

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
      std::ostringstream msg;
      msg << "Bad dimension '" << dimension() << "' in elastic material." << std::endl;
      throw std::logic_error(msg.str());
    } // switch
  
  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
    
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    // Dimensionalize coordinates for querying
    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);
    
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
	msg << ") in material '" << label() << "' using spatial database '" << _dbInitialStrain->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
    } // for

    strainVisitor.setClosure(&strainCell[0], strainCell.size(), cell, ADD_VALUES);
  } // for

  // Close databases
  _dbInitialStrain->close();

  PYLITH_METHOD_END;
} // _initializeInitialStrain

// ----------------------------------------------------------------------
// Update stateVars (for next time step).
void
pylith::materials::ElasticMaterial::_updateStateVars(PylithScalar* const stateVars,
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
