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

#include "Material.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector

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
pylith::materials::Material::Material(const int dimension,
				      const int tensorSize,
				      const Metadata& metadata) :
  _dt(0.0),
  _properties(0),
  _stateVars(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _numPropsQuadPt(0),
  _numVarsQuadPt(0),
  _dimension(dimension),
  _tensorSize(tensorSize),
  _needNewJacobian(false),
  _isJacobianSymmetric(true),
  _dbProperties(0),
  _dbInitialState(0),
  _id(0),
  _label(""),
  _metadata(metadata)
{ // constructor
  const int numProperties = metadata.numProperties();
  for (int i=0; i < numProperties; ++i) 
    _numPropsQuadPt += metadata.getProperty(i).fiberDim;
  assert(_numPropsQuadPt >= 0);

  const int numStateVars = metadata.numStateVars();
  for (int i=0; i < numStateVars; ++i)
    _numVarsQuadPt += metadata.getStateVar(i).fiberDim;
  assert(_numVarsQuadPt >= 0);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Material::deallocate(void)
{ // deallocate
  delete _normalizer; _normalizer = 0;
  delete _properties; _properties = 0;
  delete _stateVars; _stateVars = 0;

  _dbProperties = 0; // :TODO: Use shared pointer.
  _dbInitialState = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::materials::Material::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(const topology::Mesh& mesh,
					feassemble::Quadrature<topology::Mesh>* quadrature)
{ // initialize
  assert(0 != _dbProperties);
  assert(0 != quadrature);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("MaterialsFields");

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int spaceDim = quadrature->spaceDim();

  // Get cells associated with material
  DM              dmMesh = mesh.dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", _id, &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  // Create field to hold physical properties.
  delete _properties; _properties = new topology::Field<topology::Mesh>(mesh);
  _properties->label("properties");
  assert(0 != _properties);
  int fiberDim = numQuadPts * _numPropsQuadPt;
  int_array cellsTmp(cells, numCells);

  _properties->newSection(cellsTmp, fiberDim);
  _properties->allocate();
  _properties->zero();
  PetscSection propertiesSection = _properties->petscSection();
  PetscVec propertiesVec = _properties->localVector();

#if !defined(PRECOMPUTE_GEOMETRY)
  scalar_array coordinatesCell(numBasis*spaceDim);
  PetscSection coordSection;
  PetscVec coordVec;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);
#endif

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array propertiesQuery(numDBProperties);
  scalar_array propertiesCell(numQuadPts*_numPropsQuadPt);

  // Setup database for quering for physical properties
  assert(0 != _dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  // Create field to hold state variables. We create the field even
  // if there is no initial state, because this we will use this field
  // to hold the state variables.
  PetscSection stateVarsSection = PETSC_NULL;
  PetscVec stateVarsVec = PETSC_NULL;
  delete _stateVars; _stateVars = new topology::Field<topology::Mesh>(mesh);
  _stateVars->label("state variables");
  fiberDim = numQuadPts * _numVarsQuadPt;
  if (fiberDim > 0) {
    assert(0 != _stateVars);
    assert(0 != _properties);
    _stateVars->newSection(*_properties, fiberDim);
    _stateVars->allocate();
    _stateVars->zero();
    stateVarsSection = _stateVars->petscSection();
    stateVarsVec     = _stateVars->localVector();
    assert(stateVarsSection);assert(stateVarsVec);
  } // if

  // Create arrays for querying
  const int numDBStateVars = _metadata.numDBStateVars();
  scalar_array stateVarsQuery;
  scalar_array stateVarsCell;
  if (0 != _dbInitialState) {
    assert(stateVarsSection);
    assert(numDBStateVars > 0);
    assert(_numVarsQuadPt > 0);
    stateVarsQuery.resize(numDBStateVars);
    stateVarsCell.resize(numQuadPts*numDBStateVars);
    
    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(),
			       _metadata.numDBStateVars());
  } // if

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  PetscScalar *propertiesArray, *stateVarsArray;

  err = VecGetArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);
  if (stateVarsVec) {err = VecGetArray(stateVarsVec,  &stateVarsArray);CHECK_PETSC_ERROR(err);}
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(cell);
#else
    const PetscScalar *coords;
    PetscInt           coordsSize;
    err = DMPlexVecGetClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
    for(PetscInt i = 0; i < coordsSize; ++i) {coordinatesCell[i] = coords[i];}
    quadrature->computeGeometry(coordinatesCell, cell);
    err = DMPlexVecRestoreClosure(dmMesh, coordSection, coordVec, cell, &coordsSize, &coords);CHECK_PETSC_ERROR(err);
#endif

    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
				lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; 
	 iQuadPt < numQuadPts; 
	 ++iQuadPt, index+=spaceDim) {
      int err = _dbProperties->query(&propertiesQuery[0], numDBProperties,
				     &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[index+i];
	msg << ") in material " << _label << "\n"
	    << "using spatial database '" << _dbProperties->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesCell[iQuadPt*_numPropsQuadPt], 
		      propertiesQuery);
      _nondimProperties(&propertiesCell[iQuadPt*_numPropsQuadPt],
			_numPropsQuadPt);

      if (0 != _dbInitialState) {
	err = _dbInitialState->query(&stateVarsQuery[0], numDBStateVars,
				     &quadPtsGlobal[index], spaceDim, cs);
	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial state variables at \n" << "(";
	  for (int i=0; i < spaceDim; ++i)
	    msg << "  " << quadPtsGlobal[index+i];
	  msg << ") in material " << _label << "\n"
	      << "using spatial database '" << _dbInitialState->label()
	      << "'.";
	  throw std::runtime_error(msg.str());
	} // if
	_dbToStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt], 
		       stateVarsQuery);
	_nondimStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt],
			 _numVarsQuadPt);
      } // if

    } // for
    // Insert cell contribution into fields
    PetscInt dof, off, d;

    err = PetscSectionGetDof(propertiesSection, cell, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(propertiesSection, cell, &off);CHECK_PETSC_ERROR(err);
    assert(dof == numQuadPts*_numPropsQuadPt);
    for(PetscInt d = 0; d < dof; ++d) {
      propertiesArray[off+d] = propertiesCell[d];
    }
    if (0 != _dbInitialState) {
      err = PetscSectionGetDof(stateVarsSection, cell, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(stateVarsSection, cell, &off);CHECK_PETSC_ERROR(err);
      assert(dof == numQuadPts*numDBStateVars);
      for(PetscInt d = 0; d < dof; ++d) {
        stateVarsArray[off+d] = stateVarsCell[d];
      }
    }
  } // for
  err = VecRestoreArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);
  if (stateVarsVec) {err = VecRestoreArray(stateVarsVec,  &stateVarsArray);CHECK_PETSC_ERROR(err);}
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

  // Close databases
  _dbProperties->close();
  if (0 != _dbInitialState)
    _dbInitialState->close();

  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Get the properties field.
const pylith::topology::Field<pylith::topology::Mesh>*
pylith::materials::Material::propertiesField() const
{ // propertiesField
  return _properties;
} // propertiesField

// ----------------------------------------------------------------------
// Get the state variables field.
const pylith::topology::Field<pylith::topology::Mesh>*
pylith::materials::Material::stateVarsField() const
{ // stateVarsField
  return _stateVars;
} // stateVarsField

// ----------------------------------------------------------------------
// Check whether material has a field as a property.
bool
pylith::materials::Material::hasProperty(const char* name)
{ // hasProperty
  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  return (propertyIndex >= 0);
} // hasProperty

// ----------------------------------------------------------------------
// Check whether material has a field as a state variable.
bool
pylith::materials::Material::hasStateVar(const char* name)
{ // hasStateVar
  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  return (stateVarIndex >= 0);
} // hasStateVar

// ----------------------------------------------------------------------
// Get physical property or state variable field.
void
pylith::materials::Material::getField(topology::Field<topology::Mesh> *field,
				      const char* name) const
{ // getField
  // Logging of allocation is handled by getField() caller since it
  // manages the memory for the field argument.

  assert(0 != field);
  assert(0 != _properties);
  assert(0 != _stateVars);

  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);
  if (propertyIndex < 0 && stateVarIndex < 0) {
    std::ostringstream msg;
    msg << "Unknown physical property or state variable '" << name
	<< "' for material '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Get cell information
  DM              dmMesh = field->mesh().dmMesh();
  IS              cellIS;
  const PetscInt *cells;
  PetscInt        numCells;
  PetscErrorCode  err;

  assert(dmMesh);
  err = DMPlexGetStratumIS(dmMesh, "material-id", _id, &cellIS);CHECK_PETSC_ERROR(err);
  err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
  err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);

  topology::FieldBase::VectorFieldEnum fieldType = topology::FieldBase::OTHER;
      
  if (propertyIndex >= 0) { // If field is a property
    int propOffset = 0;
    assert(propertyIndex < _metadata.numProperties());
    for (int i=0; i < propertyIndex; ++i)
      propOffset += _metadata.getProperty(i).fiberDim;
    const int fiberDim = _metadata.getProperty(propertyIndex).fiberDim;

    // Get properties section
    PetscSection propertiesSection = _properties->petscSection();
    Vec          propertiesVec     = _properties->localVector();
    assert(propertiesSection);
    PetscInt totalPropsFiberDimLocal = 0;
    PetscInt totalPropsFiberDim = 0;
    if (numCells > 0) {err = PetscSectionGetDof(propertiesSection, cells[0], &totalPropsFiberDimLocal);CHECK_PETSC_ERROR(err);}
    MPI_Allreduce((void *) &totalPropsFiberDimLocal, 
                  (void *) &totalPropsFiberDim, 1, 
                  MPIU_INT, MPI_MAX, field->mesh().comm());
    assert(totalPropsFiberDim > 0);
    const int numPropsQuadPt = _numPropsQuadPt;
    const int numQuadPts = totalPropsFiberDim / numPropsQuadPt;
    assert(totalPropsFiberDim == numQuadPts * numPropsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for property field if necessary.
    PetscSection fieldSection    = field->petscSection();
    bool         useCurrentField = PETSC_FALSE;
    PetscInt     pStart, pEnd;

    err = PetscSectionGetChart(fieldSection, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
    if (pEnd < 0) {
      err = DMPlexGetHeightStratum(dmMesh, 0, &pStart, &pEnd);CHECK_PETSC_ERROR(err);
      err = PetscSectionSetChart(fieldSection, pStart, pEnd);CHECK_PETSC_ERROR(err);
    } else {
      // check fiber dimension
      PetscInt totalFiberDimCurrentLocal = 0;
      PetscInt totalFiberDimCurrent = 0;
      if (numCells > 0) {err = PetscSectionGetDof(fieldSection, cells[0], &totalFiberDimCurrentLocal);CHECK_PETSC_ERROR(err);}
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal, 
                    (void *) &totalFiberDimCurrent, 1, 
                    MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      int_array cellsTmp(cells, numCells);

      field->newSection(cellsTmp, totalFiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    assert(fieldSection);
    field->label(name);
    field->scale(1.0);
    fieldType = _metadata.getProperty(propertyIndex).fieldType;

    // Buffer for property at cell's quadrature points
    scalar_array propertiesCell(numPropsQuadPt);

    // Loop over cells
    Vec          fieldVec = field->localVector();
    PetscScalar *fieldArray, *propertiesArray;

    err = VecGetArray(fieldVec,      &fieldArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);
    for(PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];
      PetscInt       off, poff;
   
      err = PetscSectionGetOffset(fieldSection,      cell, &off);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(propertiesSection, cell, &poff);CHECK_PETSC_ERROR(err);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
        for (int i=0; i < numPropsQuadPt; ++i)
          propertiesCell[i] = propertiesArray[iQuad*numPropsQuadPt + poff+i];
        _dimProperties(&propertiesCell[0], numPropsQuadPt);
        for (int i=0; i < fiberDim; ++i)
          fieldArray[iQuad*fiberDim + off+i] = propertiesCell[propOffset+i];
      } // for
    } // for
    err = VecRestoreArray(fieldVec, &fieldArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(propertiesVec, &propertiesArray);CHECK_PETSC_ERROR(err);
  } else { // field is a state variable
    assert(stateVarIndex >= 0);
    
    int varOffset = 0;
    assert(stateVarIndex < _metadata.numStateVars());
    for (int i=0; i < stateVarIndex; ++i)
      varOffset += _metadata.getStateVar(i).fiberDim;
    const int fiberDim = _metadata.getStateVar(stateVarIndex).fiberDim;

    // Get state variables section
    PetscSection stateVarsSection = _stateVars->petscSection();
    Vec          stateVarsVec     = _stateVars->localVector();
    assert(stateVarsSection);assert(stateVarsVec);
    PetscInt totalVarsFiberDimLocal = 0;
    PetscInt totalVarsFiberDim = 0;
    if (numCells > 0) {err = PetscSectionGetDof(stateVarsSection, cells[0], &totalVarsFiberDimLocal);CHECK_PETSC_ERROR(err);}
    MPI_Allreduce((void *) &totalVarsFiberDimLocal, 
                  (void *) &totalVarsFiberDim, 1, 
                  MPIU_INT, MPI_MAX, field->mesh().comm());
    assert(totalVarsFiberDim > 0);
    const int numVarsQuadPt = _numVarsQuadPt;
    const int numQuadPts = totalVarsFiberDim / numVarsQuadPt;
    assert(totalVarsFiberDim == numQuadPts * numVarsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for state variable field if necessary.
    PetscSection fieldSection    = field->petscSection();
    PetscVec          fieldVec        = field->localVector();
    bool useCurrentField = fieldSection != PETSC_NULL;
    if (fieldSection) {
      // check fiber dimension
      PetscInt totalFiberDimCurrentLocal = 0;
      PetscInt totalFiberDimCurrent = 0;
      if (numCells > 0) {err = PetscSectionGetDof(fieldSection, cells[0], &totalFiberDimCurrentLocal);CHECK_PETSC_ERROR(err);}
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal, 
                    (void *) &totalFiberDimCurrent, 1, 
                    MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
      logger.stagePush("OutputFields");
      int_array cellsTmp(cells, numCells);

      field->newSection(cellsTmp, totalFiberDim);
      field->allocate();
      logger.stagePop();
    } // if
    assert(fieldSection);
    fieldType = _metadata.getStateVar(stateVarIndex).fieldType;
    field->label(name);
    field->scale(1.0);

    // Buffer for state variable at cell's quadrature points
    scalar_array stateVarsCell(numVarsQuadPt);
    
    // Loop over cells
    PetscScalar  *fieldArray, *stateVarsArray;

    err = VecGetArray(fieldVec,     &fieldArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(stateVarsVec, &stateVarsArray);CHECK_PETSC_ERROR(err);
    for(PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];
      PetscInt       off, soff;
      
      err = PetscSectionGetOffset(fieldSection,     cell, &off);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(stateVarsSection, cell, &soff);CHECK_PETSC_ERROR(err);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
        for (int i=0; i < numVarsQuadPt; ++i)
          stateVarsCell[i] = stateVarsArray[iQuad*numVarsQuadPt + soff+i];
        _dimStateVars(&stateVarsCell[0], numVarsQuadPt);
        for (int i=0; i < fiberDim; ++i)
          fieldArray[iQuad*fiberDim + off+i] = stateVarsCell[varOffset+i];
      } // for
    } // for
    err = VecRestoreArray(fieldVec,     &fieldArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(stateVarsVec, &stateVarsArray);CHECK_PETSC_ERROR(err);
  } // if/else
  err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
  err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);

  topology::FieldBase::VectorFieldEnum multiType = 
    topology::FieldBase::MULTI_OTHER;
  switch (fieldType)
    { // switch
    case topology::FieldBase::SCALAR:
      multiType = topology::FieldBase::MULTI_SCALAR;
      break;
    case topology::FieldBase::VECTOR:
      multiType = topology::FieldBase::MULTI_VECTOR;
      break;
    case topology::FieldBase::TENSOR:
      multiType = topology::FieldBase::MULTI_TENSOR;
      break;
    case topology::FieldBase::OTHER:
      multiType = topology::FieldBase::MULTI_OTHER;
      break;
    case topology::FieldBase::MULTI_SCALAR:
    case topology::FieldBase::MULTI_VECTOR:
    case topology::FieldBase::MULTI_TENSOR:
    case topology::FieldBase::MULTI_OTHER:
    default :
      std::cerr << "Bad vector field type '" << fieldType << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type for Material.");
    } // switch
  field->vectorFieldType(multiType);
} // getField
  
// ----------------------------------------------------------------------
// Get indices for physical property or state variable field.
void
pylith::materials::Material::_findField(int* propertyIndex,
					int* stateVarIndex,
					const char* name) const
{ // _findField
  assert(0 != propertyIndex);
  assert(0 != stateVarIndex);

  *propertyIndex = -1;
  *stateVarIndex = -1;

  const std::string nameString = name;
  const int numProperties = _metadata.numProperties();
  for (int i=0; i < numProperties; ++i)
    if (nameString == _metadata.getProperty(i).name) {
      *propertyIndex = i;
      return;
    } // if

  const int numStateVars = _metadata.numStateVars();
  for (int i=0; i < numStateVars; ++i)
    if (nameString == _metadata.getStateVar(i).name) {
      *stateVarIndex = i;
      return;
    } // if
} // _findField
  

// End of file 
