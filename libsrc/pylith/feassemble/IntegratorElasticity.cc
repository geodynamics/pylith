// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "IntegratorElasticity.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/array.hh" // USES scalar_array

#include <cstring> // USES memcpy()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorElasticity::IntegratorElasticity(void) :
  _material(0),
  _outputFields(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorElasticity::~IntegratorElasticity(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorElasticity::deallocate(void)
{ // deallocate
  Integrator<Quadrature<topology::Mesh> >::deallocate();

  delete _outputFields; _outputFields = 0;
  _material = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set material.
void
pylith::feassemble::IntegratorElasticity::material(materials::ElasticMaterial* m)
{ // material
  _material = m;
  if (0 != _material)
    _material->timeStep(_dt);  
} // material

// ----------------------------------------------------------------------
// Determine whether we need to recompute the Jacobian.
bool
pylith::feassemble::IntegratorElasticity::needNewJacobian(void)
{ // needNewJacobian
  assert(0 != _material);
  if (!_needNewJacobian)
    _needNewJacobian = _material->needNewJacobian();
  return _needNewJacobian;
} // needNewJacobian

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
void
pylith::feassemble::IntegratorElasticity::useSolnIncr(const bool flag)
{ // useSolnIncr
  Integrator<Quadrature<topology::Mesh> >::useSolnIncr(flag);
  
  assert(0 != _material);
  _material->useElasticBehavior(!flag);
} // useSolnIncr

// ----------------------------------------------------------------------
// Initialize integrator.
void
pylith::feassemble::IntegratorElasticity::initialize(const topology::Mesh& mesh)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _material);

  _initializeLogger();

  // Compute geometry for quadrature operations.
  _quadrature->initializeGeometry();
#if defined(PRECOMPUTE_GEOMETRY)
  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  _quadrature->computeGeometry(mesh, cells);
#endif

  // Initialize material.
  _material->initialize(mesh, _quadrature);
  _isJacobianSymmetric = _material->isJacobianSymmetric();

  // Allocate vectors and matrices for cell values.
  _initCellVector();
  _initCellMatrix();

  // Set up gravity field database for querying
  if (0 != _gravityField) {
    const int spaceDim = _quadrature->spaceDim();
    _gravityField->open();
    if (1 == spaceDim){
      const char* queryNames[] = { "x"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else if (2 == spaceDim){
      const char* queryNames[] = { "x", "y"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else if (3 == spaceDim){
      const char* queryNames[] = { "x", "y", "z"};
      _gravityField->queryVals(queryNames, spaceDim);
    } else {
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in IntegratorElasticity.");
    } // else
  } // if
} // initialize

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorElasticity::updateStateVars(
				      const PylithScalar t,
				      topology::SolutionFields* const fields)
{ // updateState
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != fields);

  // No need to update state vars if material doesn't have any.
  if (!_material->hasStateVars())
    return;

  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in "
			     "IntegratorElasticity::updateStateVars().");
  } // else

  // Allocate arrays for cell data.
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get fields
  scalar_array dispTCell(numBasis*spaceDim);
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMComplexGetCoordinateVec(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Retrieve geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(cell);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Restrict input fields to cell
    const PetscScalar *dispTArray;
    PetscInt           dispTSize;
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < dispTSize; ++i) {dispTCell[i] = dispTArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);

    // Get cell geometry information that depends on cell
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispTCell, numBasis, numQuadPts);

    // Update material state
    _material->updateStateVars(strainCell, cell);
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
} // updateStateVars

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorElasticity::verifyConfiguration(
					   const topology::Mesh& mesh) const
{ // verifyConfiguration
  assert(0 != _quadrature);
  assert(0 != _material);

  const int dimension = mesh.dimension();

  // check compatibility of mesh and material
  if (_material->dimension() != dimension) {
    std::ostringstream msg;
    msg << "Material '" << _material->label()
	<< "' is incompatible with mesh.\n"
	<< "Dimension of mesh: " << dimension
	<< ", dimension of material: " << _material->dimension()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if

  // check compatibility of mesh and quadrature scheme
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for material '"
	<< _material->label() << "'.\n"
	<< "Dimension of mesh: " << dimension
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->refGeometry().numCorners();

  DM dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    PetscInt       cellNumCorners;

    err = DMComplexGetConeSize(dmMesh, cell, &cellNumCorners);CHECK_PETSC_ERROR(err);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell in material '"
	  << _material->label() << "'.\n"
	  << "Cell " << cell << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::Mesh>&
pylith::feassemble::IntegratorElasticity::cellField(
					   const char* name,
					   const topology::Mesh& mesh,
					   topology::SolutionFields* fields)
{ // cellField
  assert(0 != _material);
  assert(0 != _normalizer);

  if (0 == _outputFields)
    _outputFields =
      new topology::Fields<topology::Field<topology::Mesh> >(mesh);
  
  if (0 == strcasecmp(name, "total_strain")) {

    if (_material->hasStateVar("total_strain")) {
      // total strain is available as a state variable
      assert(0 != fields);
      _allocateTensorField(mesh);
      topology::Field<topology::Mesh>& buffer = 
        _outputFields->get("buffer (tensor)");    
      _material->getField(&buffer, "total_strain");
      buffer.addDimensionOkay(true);
      return buffer;

    } else { // must calculate total strain
      assert(0 != fields);
      _allocateTensorField(mesh);
      topology::Field<topology::Mesh>& buffer = 
        _outputFields->get("buffer (tensor)");
      buffer.label("total_strain");
      buffer.scale(1.0);
      buffer.addDimensionOkay(true);
      _calcStrainStressField(&buffer, name, fields);
      return buffer;

    } // if/else
  } else if (0 == strcasecmp(name, "stress")) {

    if (_material->hasStateVar("stress")) {
      // stress is available as a state variable
      assert(0 != fields);
      _allocateTensorField(mesh);
      topology::Field<topology::Mesh>& buffer = 
        _outputFields->get("buffer (tensor)");    
      _material->getField(&buffer, "stress");
      buffer.addDimensionOkay(true);
      return buffer;

    } else { // must calculate stress from strain
      if (_material->hasStateVar("strain")) {
	// total strain is available as a state variable
	assert(0 != fields);
	_allocateTensorField(mesh);
	topology::Field<topology::Mesh>& buffer = 
	  _outputFields->get("buffer (tensor)");    
	_material->getField(&buffer, "total_strain");
	buffer.label(name);
	buffer.scale(_normalizer->pressureScale());
	buffer.addDimensionOkay(true);
	_calcStressFromStrain(&buffer);
	return buffer;

      } else { // must calculate strain 
	assert(0 != fields);
	_allocateTensorField(mesh);
	topology::Field<topology::Mesh>& buffer = 
	  _outputFields->get("buffer (tensor)");
	buffer.label("stress");
	buffer.scale(_normalizer->pressureScale());
	buffer.addDimensionOkay(true);
	_calcStrainStressField(&buffer, name, fields);
	return buffer;
	
      } // else

    } // else
  } else {
    if (!_outputFields->hasField("buffer (other)"))
      _outputFields->add("buffer (other)", "buffer");
    topology::Field<topology::Mesh>& buffer =
      _outputFields->get("buffer (other)");
    _material->getField(&buffer, name);
    buffer.addDimensionOkay(true);
    return buffer;
    
  } // if/else

  
  // Return tensor section to satisfy member function definition. Code
  // should never get here.
  throw std::logic_error("Internal error.");
  topology::Field<topology::Mesh>& buffer = 
    _outputFields->get("buffer (tensor)");    

  return buffer;
} // cellField

// ----------------------------------------------------------------------
// Get output fields.
const pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >*
pylith::feassemble::IntegratorElasticity::outputFields(void) const
{ // outputFields
  return _outputFields;
} // outputFields

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::feassemble::IntegratorElasticity::_initializeLogger(void)
{ // initializeLogger
  delete _logger; _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("ElasticityIntegrator");
  _logger->initialize();
  _logger->registerEvent("ElIR setup");
  _logger->registerEvent("ElIR geometry");
  _logger->registerEvent("ElIR compute");
  _logger->registerEvent("ElIR restrict");
  _logger->registerEvent("ElIR stateVars");
  _logger->registerEvent("ElIR stress");
  _logger->registerEvent("ElIR update");
 
  _logger->registerEvent("ElIJ setup");
  _logger->registerEvent("ElIJ geometry");
  _logger->registerEvent("ElIJ compute");
  _logger->registerEvent("ElIJ restrict");
  _logger->registerEvent("ElIJ stateVars");
  _logger->registerEvent("ElIJ update");
} // initializeLogger

// ----------------------------------------------------------------------
// Allocate buffer for tensor field at quadrature points.
void
pylith::feassemble::IntegratorElasticity::_allocateTensorField(
						 const topology::Mesh& mesh)
{ // _allocateTensorField
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _outputFields);

  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  int_array cellsTmp(cells, numCells);

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  
  if (!_outputFields->hasField("buffer (tensor)")) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    logger.stagePush("OutputFields");
    _outputFields->add("buffer (tensor)", "buffer");
    topology::Field<topology::Mesh>& buffer =
      _outputFields->get("buffer (tensor)");
    buffer.newSection(cellsTmp, numQuadPts*tensorSize);
    buffer.allocate();
    buffer.vectorFieldType(topology::FieldBase::MULTI_TENSOR);
    buffer.addDimensionOkay(true);
    logger.stagePop();
  } // if
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
} // _allocateTensorField

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcStrainStressField(
				 topology::Field<topology::Mesh>* field,
				 const char* name,
				 topology::SolutionFields* const fields)
{ // _calcStrainStressField
  assert(0 != field);
  assert(0 != _quadrature);
  assert(0 != _material);

  const bool calcStress = (0 == strcasecmp(name, "stress")) ? true : false;
    
  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in IntegratorElasticity.");
  } // else
  
  // Allocate arrays for cell data.
  scalar_array dispCell(numBasis*spaceDim);
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  scalar_array stressCell(numQuadPts*tensorSize);
  stressCell = 0.0;

  // Get cell information
  DM              dmMesh = fields->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get field
  scalar_array dispTCell(numBasis*spaceDim);
  topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
  PetscSection dispTSection = dispT.petscSection();
  Vec          dispTVec     = dispT.localVector();
  assert(dispTSection);assert(dispTVec);
    
#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMComplexGetCoordinateVec(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  PetscSection fieldSection = field->petscSection();
  assert(fieldSection);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Retrieve geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMComplexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    _quadrature->computeGeometry(coordinatesCell, cell);
    err = DMComplexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    // Restrict input fields to cell
    const PetscScalar *dispTArray;
    PetscInt           dispTSize;
    err = DMComplexVecGetClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < dispTSize; ++i) {dispTCell[i] = dispTArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, dispTSection, dispTVec, cell, &dispTSize, &dispTArray);CHECK_PETSC_ERROR(err);

    // Get cell geometry information that depends on cell
    const scalar_array& basisDeriv = _quadrature->basisDeriv();
    
    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispTCell, numBasis, numQuadPts);
    
    if (!calcStress) {
      err = DMComplexVecSetClosure(dmMesh, fieldSection, PETSC_NULL, cell, &strainCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
    } else {
      _material->retrievePropsAndVars(cell);
      stressCell = _material->calcStress(strainCell);
      err = DMComplexVecSetClosure(dmMesh, fieldSection, PETSC_NULL, cell, &stressCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
    } // else
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
} // _calcStrainStressField

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcStressFromStrain(
				   topology::Field<topology::Mesh>* field)
{ // _calcStressFromStrain
  assert(0 != field);
  assert(0 != _quadrature);
  assert(0 != _material);

  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int tensorSize = _material->tensorSize();
  
  // Allocate arrays for cell data.
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  scalar_array stressCell(numQuadPts*tensorSize);
  stressCell = 0.0;

  // Get cell information
  DM              dmMesh = field->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMComplexGetStratumIS(dmMesh, "material-id", _material->id(), &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  // Get field
  PetscSection fieldSection = field->petscSection();
  Vec          fieldVec     = field->localVector();
  assert(fieldSection);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt     cell = cells[c];
    PetscInt           fieldSize;
    const PetscScalar *fieldArray;

    err = DMComplexVecGetClosure(dmMesh, fieldSection, fieldVec, cell, &fieldSize, &fieldArray);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < fieldSize; ++i) {strainCell[i] = fieldArray[i];}
    err = DMComplexVecRestoreClosure(dmMesh, fieldSection, fieldVec, cell, &fieldSize, &fieldArray);CHECK_PETSC_ERROR(err);

    _material->retrievePropsAndVars(cell);
    stressCell = _material->calcStress(strainCell);
    err = DMComplexVecSetClosure(dmMesh, fieldSection, PETSC_NULL, cell, &stressCell[0], ADD_VALUES);CHECK_PETSC_ERROR(err);
  } // for
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
} // _calcStressFromStrain

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 1-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual1D(
				     const scalar_array& stress)
{ // _elasticityResidual1D
  const int spaceDim = 1;
  const int cellDim = 1;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();

  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar s11 = stress[iQuad];
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const PylithScalar N1 = wt*basisDeriv[iQuad*numBasis+iBasis  ];
      _cellVector[iBasis*spaceDim  ] -= N1*s11;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*5));
} // _elasticityResidual1D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 2-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual2D(
				     const scalar_array& stress)
{ // _elasticityResidual2D
  const int cellDim = 2;
  const int spaceDim = 2;
  const int stressSize = 3;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQs = iQuad*stressSize;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar s11 = stress[iQs  ];
    const PylithScalar s22 = stress[iQs+1];
    const PylithScalar s12 = stress[iQs+2];
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const int iBlock = iBasis*spaceDim;
      const PylithScalar N1 = wt*basisDeriv[iQ+iBlock  ];
      const PylithScalar N2 = wt*basisDeriv[iQ+iBlock+1];

      _cellVector[iBlock  ] -= N1*s11 + N2*s12;
      _cellVector[iBlock+1] -= N1*s12 + N2*s22;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(8+2+9)));
} // _elasticityResidual2D

// ----------------------------------------------------------------------
// Integrate elasticity term in residual for 3-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityResidual3D(
				     const scalar_array& stress)
{ // _elasticityResidual3D
  const int spaceDim = 3;
  const int cellDim = 3;
  const int stressSize = 6;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const int iQs = iQuad * stressSize;
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar s11 = stress[iQs  ];
    const PylithScalar s22 = stress[iQs+1];
    const PylithScalar s33 = stress[iQs+2];
    const PylithScalar s12 = stress[iQs+3];
    const PylithScalar s23 = stress[iQs+4];
    const PylithScalar s13 = stress[iQs+5];
    
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
        iBasis < numBasis;
        ++iBasis) {
      const int iBlock = iBasis*spaceDim;
      const PylithScalar N1 = wt*basisDeriv[iQ+iBlock+0];
      const PylithScalar N2 = wt*basisDeriv[iQ+iBlock+1];
      const PylithScalar N3 = wt*basisDeriv[iQ+iBlock+2];

      _cellVector[iBlock  ] -= N1*s11 + N2*s12 + N3*s13;
      _cellVector[iBlock+1] -= N1*s12 + N2*s22 + N3*s23;
      _cellVector[iBlock+2] -= N1*s13 + N2*s23 + N3*s33;
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+12)));
} // _elasticityResidual3D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 1-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian1D(
			       const scalar_array& elasticConsts)
{ // _elasticityJacobian1D
  const int cellDim = 1;
  const int spaceDim = 1;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);
  
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    const PylithScalar C1111 = elasticConsts[iQuad];
    for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
      const PylithScalar valI = wt*basisDeriv[iQ+iBasis]*C1111;
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const PylithScalar valIJ = valI * basisDeriv[iQ+jBasis];
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	_cellMatrix[iBlock+jBlock] += valIJ;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*3)));
} // _elasticityJacobian1D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 2-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian2D(
			       const scalar_array& elasticConsts)
{ // _elasticityJacobian2D
  const int spaceDim = 2;
  const int cellDim = 2;
  const int numConsts = 9;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);

  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const PylithScalar C1111 = elasticConsts[iQuad*numConsts+0];
    const PylithScalar C1122 = elasticConsts[iQuad*numConsts+1];
    const PylithScalar C1112 = elasticConsts[iQuad*numConsts+2]/2.0;
    const PylithScalar C2211 = elasticConsts[iQuad*numConsts+3];
    const PylithScalar C2222 = elasticConsts[iQuad*numConsts+4];
    const PylithScalar C2212 = elasticConsts[iQuad*numConsts+5]/2.0;
    const PylithScalar C1211 = elasticConsts[iQuad*numConsts+6];
    const PylithScalar C1222 = elasticConsts[iQuad*numConsts+7];
    const PylithScalar C1212 = elasticConsts[iQuad*numConsts+8]/2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const PylithScalar Ni1 = wt*basisDeriv[iQ+iBasis*spaceDim  ];
      const PylithScalar Ni2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const int iBlock = (iBasis*spaceDim  ) * (numBasis*spaceDim);
      const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const PylithScalar Nj1 = basisDeriv[iQ+jBasis*spaceDim  ];
	const PylithScalar Nj2 = basisDeriv[iQ+jBasis*spaceDim+1];
	const PylithScalar ki0j0 = 
	  C1111 * Ni1 * Nj1 + C1211 * Ni2 * Nj1 +
	  C1112 * Ni1 * Nj2 + C1212 * Ni2 * Nj2;
	const PylithScalar ki0j1 =
	  C1122 * Ni1 * Nj2 + C1222 * Ni2 * Nj2 +
	  C1112 * Ni1 * Nj1 + C1212 * Ni2 * Nj1;
	const PylithScalar ki1j0 =
	  C2211 * Ni2 * Nj1 + C1211 * Ni1 * Nj1 +
	  C2212 * Ni2 * Nj2 + C1212 * Ni1 * Nj2;
	const PylithScalar ki1j1 =
	  C2222 * Ni2 * Nj2 + C1222 * Ni1 * Nj2 +
	  C2212 * Ni2 * Nj1 + C1212 * Ni1 * Nj1;
	const int jBlock = (jBasis*spaceDim  );
	const int jBlock1 = (jBasis*spaceDim+1);
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(2+numBasis*(3*11+4))));
} // _elasticityJacobian2D

// ----------------------------------------------------------------------
// Integrate elasticity term in Jacobian for 3-D cells.
void
pylith::feassemble::IntegratorElasticity::_elasticityJacobian3D(
			       const scalar_array& elasticConsts)
{ // _elasticityJacobian3D
  const int spaceDim = 3;
  const int cellDim = 3;
  const int numConsts = 36;

  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const scalar_array& quadWts = _quadrature->quadWts();
  const scalar_array& jacobianDet = _quadrature->jacobianDet();
  const scalar_array& basisDeriv = _quadrature->basisDeriv();
  
  assert(_quadrature->spaceDim() == spaceDim);
  assert(_quadrature->cellDim() == cellDim);
  assert(quadWts.size() == numQuadPts);

  // Compute Jacobian for consistent tangent matrix
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    const PylithScalar wt = quadWts[iQuad] * jacobianDet[iQuad];
    // tau_ij = C_ijkl * e_kl
    //        = C_ijlk * 0.5 (u_k,l + u_l,k)
    //        = 0.5 * C_ijkl * (u_k,l + u_l,k)
    // divide C_ijkl by 2 if k != l
    const PylithScalar C1111 = elasticConsts[iQuad*numConsts+ 0];
    const PylithScalar C1122 = elasticConsts[iQuad*numConsts+ 1];
    const PylithScalar C1133 = elasticConsts[iQuad*numConsts+ 2];
    const PylithScalar C1112 = elasticConsts[iQuad*numConsts+ 3] / 2.0;
    const PylithScalar C1123 = elasticConsts[iQuad*numConsts+ 4] / 2.0;
    const PylithScalar C1113 = elasticConsts[iQuad*numConsts+ 5] / 2.0;
    const PylithScalar C2211 = elasticConsts[iQuad*numConsts+ 6];
    const PylithScalar C2222 = elasticConsts[iQuad*numConsts+ 7];
    const PylithScalar C2233 = elasticConsts[iQuad*numConsts+ 8];
    const PylithScalar C2212 = elasticConsts[iQuad*numConsts+ 9] / 2.0;
    const PylithScalar C2223 = elasticConsts[iQuad*numConsts+10] / 2.0;
    const PylithScalar C2213 = elasticConsts[iQuad*numConsts+11] / 2.0;
    const PylithScalar C3311 = elasticConsts[iQuad*numConsts+12];
    const PylithScalar C3322 = elasticConsts[iQuad*numConsts+13];
    const PylithScalar C3333 = elasticConsts[iQuad*numConsts+14];
    const PylithScalar C3312 = elasticConsts[iQuad*numConsts+15] / 2.0;
    const PylithScalar C3323 = elasticConsts[iQuad*numConsts+16] / 2.0;
    const PylithScalar C3313 = elasticConsts[iQuad*numConsts+17] / 2.0;
    const PylithScalar C1211 = elasticConsts[iQuad*numConsts+18];
    const PylithScalar C1222 = elasticConsts[iQuad*numConsts+19];
    const PylithScalar C1233 = elasticConsts[iQuad*numConsts+20];
    const PylithScalar C1212 = elasticConsts[iQuad*numConsts+21] / 2.0;
    const PylithScalar C1223 = elasticConsts[iQuad*numConsts+22] / 2.0;
    const PylithScalar C1213 = elasticConsts[iQuad*numConsts+23] / 2.0;
    const PylithScalar C2311 = elasticConsts[iQuad*numConsts+24];
    const PylithScalar C2322 = elasticConsts[iQuad*numConsts+25];
    const PylithScalar C2333 = elasticConsts[iQuad*numConsts+26];
    const PylithScalar C2312 = elasticConsts[iQuad*numConsts+27] / 2.0;
    const PylithScalar C2323 = elasticConsts[iQuad*numConsts+28] / 2.0;
    const PylithScalar C2313 = elasticConsts[iQuad*numConsts+29] / 2.0;
    const PylithScalar C1311 = elasticConsts[iQuad*numConsts+30];
    const PylithScalar C1322 = elasticConsts[iQuad*numConsts+31];
    const PylithScalar C1333 = elasticConsts[iQuad*numConsts+32];
    const PylithScalar C1312 = elasticConsts[iQuad*numConsts+33] / 2.0;
    const PylithScalar C1323 = elasticConsts[iQuad*numConsts+34] / 2.0;
    const PylithScalar C1313 = elasticConsts[iQuad*numConsts+35] / 2.0;
    for (int iBasis=0, iQ=iQuad*numBasis*spaceDim;
	 iBasis < numBasis;
	 ++iBasis) {
      const PylithScalar Ni1 = wt*basisDeriv[iQ+iBasis*spaceDim+0];
      const PylithScalar Ni2 = wt*basisDeriv[iQ+iBasis*spaceDim+1];
      const PylithScalar Ni3 = wt*basisDeriv[iQ+iBasis*spaceDim+2];
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const PylithScalar Nj1 = basisDeriv[iQ+jBasis*spaceDim+0];
	const PylithScalar Nj2 = basisDeriv[iQ+jBasis*spaceDim+1];
	const PylithScalar Nj3 = basisDeriv[iQ+jBasis*spaceDim+2];
	const PylithScalar ki0j0 = 
	  C1111 * Ni1 * Nj1 + C1211 * Ni2 * Nj1 + C1311 * Ni3 * Nj1 +
	  C1112 * Ni1 * Nj2 + C1212 * Ni2 * Nj2 + C1312 * Ni3 * Nj2 +
	  C1113 * Ni1 * Nj3 + C1213 * Ni2 * Nj3 + C1313 * Ni3 * Nj3;
	const PylithScalar ki0j1 =
	  C1122 * Ni1 * Nj2 + C1222 * Ni2 * Nj2 + C1322 * Ni3 * Nj2 +
	  C1112 * Ni1 * Nj1 + C1212 * Ni2 * Nj1 + C1312 * Ni3 * Nj1 +
	  C1123 * Ni1 * Nj3 + C1223 * Ni2 * Nj3 + C1323 * Ni3 * Nj3;
	const PylithScalar ki0j2 =
	  C1133 * Ni1 * Nj3 + C1233 * Ni2 * Nj3 + C1333 * Ni3 * Nj3 +
	  C1123 * Ni1 * Nj2 + C1223 * Ni2 * Nj2 + C1323 * Ni3 * Nj2 +
	  C1113 * Ni1 * Nj1 + C1213 * Ni2 * Nj1 + C1313 * Ni3 * Nj1;
	const PylithScalar ki1j0 =
	  C2211 * Ni2 * Nj1 + C1211 * Ni1 * Nj1 + C2311 * Ni3 * Nj1 +
	  C2212 * Ni2 * Nj2 + C1212 * Ni1 * Nj2 + C2312 * Ni3 * Nj2 +
	  C2213 * Ni2 * Nj3 + C1213 * Ni1 * Nj3 + C2313 * Ni3 * Nj3;
	const PylithScalar ki1j1 =
	  C2222 * Ni2 * Nj2 + C1222 * Ni1 * Nj2 + C2322 * Ni3 * Nj2 +
	  C2212 * Ni2 * Nj1 + C1212 * Ni1 * Nj1 + C2312 * Ni3 * Nj1 +
	  C2223 * Ni2 * Nj3 + C1223 * Ni1 * Nj3 + C2323 * Ni3 * Nj3;
	const PylithScalar ki1j2 =
	  C2233 * Ni2 * Nj3 + C1233 * Ni1 * Nj3 + C2333 * Ni3 * Nj3 +
	  C2223 * Ni2 * Nj2 + C1223 * Ni1 * Nj2 + C2323 * Ni3 * Nj2 +
	  C2213 * Ni2 * Nj1 + C1213 * Ni1 * Nj1 + C2313 * Ni3 * Nj1;
	const PylithScalar ki2j0 =
	  C3311 * Ni3 * Nj1 + C2311 * Ni2 * Nj1 + C1311 * Ni1 * Nj1 +
	  C3312 * Ni3 * Nj2 + C2312 * Ni2 * Nj2 + C1312 * Ni1 * Nj2 +
	  C3313 * Ni3 * Nj3 + C2313 * Ni2 * Nj3 + C1313 * Ni1 * Nj3; 
	const PylithScalar ki2j1 =
	  C3322 * Ni3 * Nj2 + C2322 * Ni2 * Nj2 + C1322 * Ni1 * Nj2 +
	  C3312 * Ni3 * Nj1 + C2312 * Ni2 * Nj1 + C1312 * Ni1 * Nj1 +
	  C3323 * Ni3 * Nj3 + C2323 * Ni2 * Nj3 + C1323 * Ni1 * Nj3; 
	const PylithScalar ki2j2 =
	  C3333 * Ni3 * Nj3 + C2333 * Ni2 * Nj3 + C1333 * Ni1 * Nj3 +
	  C3323 * Ni3 * Nj2 + C2323 * Ni2 * Nj2 + C1323 * Ni1 * Nj2 +
	  C3313 * Ni3 * Nj1 + C2313 * Ni2 * Nj1 + C1313 * Ni1 * Nj1;
	const int iBlock = iBasis*spaceDim * (numBasis*spaceDim);
	const int iBlock1 = (iBasis*spaceDim+1) * (numBasis*spaceDim);
	const int iBlock2 = (iBasis*spaceDim+2) * (numBasis*spaceDim);
	const int jBlock = jBasis*spaceDim;
	const int jBlock1 = jBasis*spaceDim+1;
	const int jBlock2 = jBasis*spaceDim+2;
	_cellMatrix[iBlock +jBlock ] += ki0j0;
	_cellMatrix[iBlock +jBlock1] += ki0j1;
	_cellMatrix[iBlock +jBlock2] += ki0j2;
	_cellMatrix[iBlock1+jBlock ] += ki1j0;
	_cellMatrix[iBlock1+jBlock1] += ki1j1;
	_cellMatrix[iBlock1+jBlock2] += ki1j2;
	_cellMatrix[iBlock2+jBlock ] += ki2j0;
	_cellMatrix[iBlock2+jBlock1] += ki2j1;
	_cellMatrix[iBlock2+jBlock2] += ki2j2;
      } // for
    } // for
  } // for
  PetscLogFlops(numQuadPts*(1+numBasis*(3+numBasis*(6*26+9))));
} // _elasticityJacobian3D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D(
					    scalar_array* strain,
					    const scalar_array& basisDeriv,
					    const scalar_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain1D
  assert(0 != strain);

  const int dim = 1;
  
  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      (*strain)[iQuad] += basisDeriv[iQuad*numBasis+iBasis] * disp[iBasis];
  } // for
} // calcTotalStrain1D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D(
					    scalar_array* strain,
					    const scalar_array& basisDeriv,
					    const scalar_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain2D
  assert(0 != strain);
  
  const int dim = 2;
  const int strainSize = 3;

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad*strainSize+0] += 
	basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim  ];
      (*strain)[iQuad*strainSize+1] += 
	basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+1];
      (*strain)[iQuad*strainSize+2] += 
	0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
	       basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
    } // for
} // calcTotalStrain2D

// ----------------------------------------------------------------------
void
pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D(
					    scalar_array* strain,
					    const scalar_array& basisDeriv,
					    const scalar_array& disp,
					    const int numBasis,
					    const int numQuadPts)
{ // calcTotalStrain3D
  assert(0 != strain);

  const int dim = 3;
  const int strainSize = 6;

  assert(basisDeriv.size() == numQuadPts*numBasis*dim);
  assert(disp.size() == numBasis*dim);

  (*strain) = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    for (int iBasis=0, iQ=iQuad*numBasis*dim; iBasis < numBasis; ++iBasis) {
      (*strain)[iQuad*strainSize  ] += basisDeriv[iQ+iBasis*dim]
          * disp[iBasis*dim];
      (*strain)[iQuad*strainSize+1] += basisDeriv[iQ+iBasis*dim+1]
          * disp[iBasis*dim+1];
      (*strain)[iQuad*strainSize+2] += basisDeriv[iQ+iBasis*dim+2]
          * disp[iBasis*dim+2];
      (*strain)[iQuad*strainSize+3] +=
          0.5 * (basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim  ] +
                 basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+1]);
      (*strain)[iQuad*strainSize+4] +=
          0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim+1] +
                 basisDeriv[iQ+iBasis*dim+1] * disp[iBasis*dim+2]);
      (*strain)[iQuad*strainSize+5] +=
          0.5 * (basisDeriv[iQ+iBasis*dim+2] * disp[iBasis*dim  ] +
                 basisDeriv[iQ+iBasis*dim  ] * disp[iBasis*dim+2]);
    } // for
} // calcTotalStrain3D


// End of file 
