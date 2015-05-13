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

#include "Material.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int dimension,
				      const int tensorSize,
				      const Metadata& metadata) :
  _dt(0.0),
  _properties(0),
  _stateVars(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _materialIS(0),
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
  PYLITH_METHOD_BEGIN;

  delete _normalizer; _normalizer = 0;
  delete _materialIS; _materialIS = 0;
  delete _properties; _properties = 0;
  delete _stateVars; _stateVars = 0;

  _dbProperties = 0; // :TODO: Use shared pointer.
  _dbInitialState = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::materials::Material::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  PYLITH_METHOD_BEGIN;

  if (!_normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;

  PYLITH_METHOD_END;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(const topology::Mesh& mesh,
					feassemble::Quadrature* quadrature)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_dbProperties);
  assert(quadrature);

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int spaceDim = quadrature->spaceDim();

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);
  const PetscInt numCells = _materialIS->size();
  const PetscInt* cells = _materialIS->points();

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);

  // Create field to hold physical properties.
  delete _properties; _properties = new topology::Field(mesh);assert(_properties);
  _properties->label("properties");
  const int propsFiberDim = numQuadPts * _numPropsQuadPt;
  int_array cellsTmp(cells, numCells);

  _properties->newSection(cellsTmp, propsFiberDim);
  _properties->allocate();
  _properties->zeroAll();
  topology::VecVisitorMesh propertiesVisitor(*_properties);
  PetscScalar* propertiesArray = propertiesVisitor.localArray();

  scalar_array coordsCell(numBasis*spaceDim); // :KULDGE: Update numBasis to numCorners after implementing higher order
  topology::CoordsVisitor coordsVisitor(dmMesh);

  // Optimize coordinate retrieval in closure  
  topology::CoordsVisitor::optimizeClosure(dmMesh);

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  scalar_array quadPtsGlobal(numQuadPts*spaceDim);
  scalar_array propertiesQuery(numDBProperties);
  scalar_array propertiesCell(propsFiberDim);

  // Setup database for quering for physical properties
  assert(_dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  // Create field to hold state variables. We create the field even
  // if there is no initial state, because this we will use this field
  // to hold the state variables.
  delete _stateVars; _stateVars = new topology::Field(mesh);assert(_stateVars);
  _stateVars->label("state variables");
  const int stateVarsFiberDim = numQuadPts * _numVarsQuadPt;
  topology::VecVisitorMesh* stateVarsVisitor = 0;
  PetscScalar* stateVarsArray = NULL;
  if (stateVarsFiberDim > 0) {
    assert(_stateVars);
    assert(_properties);
    _stateVars->newSection(*_properties, stateVarsFiberDim);
    _stateVars->allocate();
    _stateVars->zeroAll();
    stateVarsVisitor = new topology::VecVisitorMesh(*_stateVars);
    stateVarsArray = stateVarsVisitor->localArray();
  } // if


  // Create arrays for querying
  const int numDBStateVars = _metadata.numDBStateVars();
  scalar_array stateVarsQuery;
  scalar_array stateVarsCell;
  if (_dbInitialState) {
    assert(numDBStateVars > 0);
    assert(_numVarsQuadPt > 0);
    stateVarsQuery.resize(numDBStateVars);
    stateVarsCell.resize(stateVarsFiberDim);
    
    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(), _metadata.numDBStateVars());
  } // if

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];

    // Compute geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);

    const scalar_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0; iQuadPt < numQuadPts; ++iQuadPt, index+=spaceDim) {
      int err = _dbProperties->query(&propertiesQuery[0], numDBProperties, &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at " << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[index+i];
	msg << ") in material '" << _label << "' using spatial database '" << _dbProperties->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesCell[iQuadPt*_numPropsQuadPt], propertiesQuery);
      _nondimProperties(&propertiesCell[iQuadPt*_numPropsQuadPt], _numPropsQuadPt);

      if (_dbInitialState) {
	err = _dbInitialState->query(&stateVarsQuery[0], numDBStateVars, &quadPtsGlobal[index], spaceDim, cs);
	if (err) {
	  std::ostringstream msg;
	  msg << "Could not find initial state variables at \n" << "(";
	  for (int i=0; i < spaceDim; ++i)
	    msg << "  " << quadPtsGlobal[index+i];
	  msg << ") in material '" << _label << "' using spatial database '" << _dbInitialState->label() << "'.";
	  throw std::runtime_error(msg.str());
	} // if
	_dbToStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt], stateVarsQuery);
	_nondimStateVars(&stateVarsCell[iQuadPt*_numVarsQuadPt], _numVarsQuadPt);
      } // if

    } // for
    // Insert cell contribution into fields
    const PetscInt off = propertiesVisitor.sectionOffset(cell);
    assert(propsFiberDim == propertiesVisitor.sectionDof(cell));
    for(PetscInt d = 0; d < propsFiberDim; ++d) {
      propertiesArray[off+d] = propertiesCell[d];
    } // for
    if (_dbInitialState) {
      assert(stateVarsVisitor);
      assert(stateVarsArray);
      const PetscInt off = stateVarsVisitor->sectionOffset(cell);
      assert(stateVarsFiberDim == stateVarsVisitor->sectionDof(cell));
      for(PetscInt d = 0; d < stateVarsFiberDim; ++d) {
        stateVarsArray[off+d] = stateVarsCell[d];
      } // for
    } // if
  } // for
  delete stateVarsVisitor; stateVarsVisitor = 0;

  // Close databases
  _dbProperties->close();
  if (_dbInitialState)
    _dbInitialState->close();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get the properties field.
const pylith::topology::Field*
pylith::materials::Material::propertiesField() const
{ // propertiesField
  return _properties;
} // propertiesField

// ----------------------------------------------------------------------
// Get the state variables field.
const pylith::topology::Field*
pylith::materials::Material::stateVarsField() const
{ // stateVarsField
  return _stateVars;
} // stateVarsField

// ----------------------------------------------------------------------
// Check whether material has a field as a property.
bool
pylith::materials::Material::hasProperty(const char* name)
{ // hasProperty
  PYLITH_METHOD_BEGIN;

  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);

  PYLITH_METHOD_RETURN(propertyIndex >= 0);
} // hasProperty

// ----------------------------------------------------------------------
// Check whether material has a field as a state variable.
bool
pylith::materials::Material::hasStateVar(const char* name)
{ // hasStateVar
  PYLITH_METHOD_BEGIN;

  int propertyIndex = -1;
  int stateVarIndex = -1;
  _findField(&propertyIndex, &stateVarIndex, name);

  PYLITH_METHOD_RETURN(stateVarIndex >= 0);
} // hasStateVar

// ----------------------------------------------------------------------
// Get physical property or state variable field.
void
pylith::materials::Material::getField(topology::Field *field,
				      const char* name) const
{ // getField
  PYLITH_METHOD_BEGIN;

  assert(field);
  assert(_properties);
  assert(_stateVars);

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
  PetscDM dmMesh = field->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt numCells = _materialIS->size();
  const PetscInt* cells = _materialIS->points();

  topology::FieldBase::VectorFieldEnum fieldType = topology::FieldBase::OTHER;
      
  if (propertyIndex >= 0) { // If field is a property
    int propOffset = 0;
    assert(propertyIndex < _metadata.numProperties());
    for (int i=0; i < propertyIndex; ++i)
      propOffset += _metadata.getProperty(i).fiberDim;
    const int fiberDim = _metadata.getProperty(propertyIndex).fiberDim;
    topology::VecVisitorMesh propertiesVisitor(*_properties);
    PetscScalar* propertiesArray = propertiesVisitor.localArray();

    // Get properties section
    PetscInt totalPropsFiberDimLocal = 0;
    PetscInt totalPropsFiberDim = 0;
    if (numCells > 0) {
      totalPropsFiberDimLocal = propertiesVisitor.sectionDof(cells[0]);
    } // if
    MPI_Allreduce((void *) &totalPropsFiberDimLocal, (void *) &totalPropsFiberDim, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
    assert(totalPropsFiberDim > 0);
    const int numPropsQuadPt = _numPropsQuadPt;
    const int numQuadPts = totalPropsFiberDim / numPropsQuadPt;
    assert(totalPropsFiberDim == numQuadPts * numPropsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for property field if necessary.
    bool useCurrentField = false;
    if (field->hasSection()) {
      // check fiber dimension
      PetscInt totalFiberDimCurrentLocal = 0;
      PetscInt totalFiberDimCurrent = 0;
      if (numCells > 0) {
	PetscSection fieldSection = field->localSection();
	PetscErrorCode err = PetscSectionGetDof(fieldSection, cells[0], &totalFiberDimCurrentLocal);PYLITH_CHECK_ERROR(err);
      } // if
      MPI_Allreduce((void *) &totalFiberDimCurrentLocal, (void *) &totalFiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, numCells, totalFiberDim);
      field->allocate();
    } // if
    field->label(name);
    field->scale(1.0);
    fieldType = _metadata.getProperty(propertyIndex).fieldType;
    topology::VecVisitorMesh fieldVisitor(*field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    // Buffer for property at cell's quadrature points
    scalar_array propertiesCell(numPropsQuadPt);

    // Loop over cells
    for(PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];

      const PetscInt poff = propertiesVisitor.sectionOffset(cell);
      const PetscInt foff = fieldVisitor.sectionOffset(cell);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
        for (int i=0; i < numPropsQuadPt; ++i)
          propertiesCell[i] = propertiesArray[iQuad*numPropsQuadPt + poff+i];
        _dimProperties(&propertiesCell[0], numPropsQuadPt);
        for (int i=0; i < fiberDim; ++i)
          fieldArray[iQuad*fiberDim + foff+i] = propertiesCell[propOffset+i];
      } // for
    } // for
  } else { // field is a state variable
    assert(stateVarIndex >= 0);
    
    int varOffset = 0;
    assert(stateVarIndex < _metadata.numStateVars());
    for (int i=0; i < stateVarIndex; ++i)
      varOffset += _metadata.getStateVar(i).fiberDim;
    const int fiberDim = _metadata.getStateVar(stateVarIndex).fiberDim;

    // Get state variables
    topology::VecVisitorMesh stateVarsVisitor(*_stateVars);
    PetscScalar* stateVarsArray = stateVarsVisitor.localArray();

    PetscInt totalVarsFiberDimLocal = 0;
    PetscInt totalVarsFiberDim = 0;
    if (numCells > 0) {
      totalVarsFiberDimLocal = stateVarsVisitor.sectionDof(cells[0]);
    } // if
    MPI_Allreduce((void*) &totalVarsFiberDimLocal, (void*) &totalVarsFiberDim, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
    assert(totalVarsFiberDim > 0);
    const int numVarsQuadPt = _numVarsQuadPt;
    const int numQuadPts = totalVarsFiberDim / numVarsQuadPt;
    assert(totalVarsFiberDim == numQuadPts * numVarsQuadPt);
    const int totalFiberDim = numQuadPts * fiberDim;

    // Allocate buffer for state variable field if necessary.
    bool useCurrentField = false;
    if (field->hasSection()) {
      // check fiber dimension
      PetscInt totalFiberDimCurrentLocal = 0;
      PetscInt totalFiberDimCurrent = 0;
      if (numCells > 0) {
	PetscSection fieldSection = field->localSection();
	PetscErrorCode err = PetscSectionGetDof(fieldSection, cells[0], &totalFiberDimCurrentLocal);PYLITH_CHECK_ERROR(err);
      } // if
      MPI_Allreduce((void*) &totalFiberDimCurrentLocal, (void*) &totalFiberDimCurrent, 1, MPIU_INT, MPI_MAX, field->mesh().comm());
      assert(totalFiberDimCurrent > 0);
      useCurrentField = totalFiberDim == totalFiberDimCurrent;
    } // if
    if (!useCurrentField) {
      field->newSection(cells, numCells, totalFiberDim);
      field->allocate();
    } // if
    fieldType = _metadata.getStateVar(stateVarIndex).fieldType;
    field->label(name);
    field->scale(1.0);
    topology::VecVisitorMesh fieldVisitor(*field);
    PetscScalar* fieldArray = fieldVisitor.localArray();

    // Buffer for state variable at cell's quadrature points
    scalar_array stateVarsCell(numVarsQuadPt);
    
    // Loop over cells
    for(PetscInt c = 0; c < numCells; ++c) {
      const PetscInt cell = cells[c];

      const PetscInt foff = fieldVisitor.sectionOffset(cell);
      const PetscInt soff = stateVarsVisitor.sectionOffset(cell);
      for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
        for (int i=0; i < numVarsQuadPt; ++i) {
          stateVarsCell[i] = stateVarsArray[iQuad*numVarsQuadPt + soff+i];
	} // for
	_dimStateVars(&stateVarsCell[0], numVarsQuadPt);
        for (int i=0; i < fiberDim; ++i)
          fieldArray[iQuad*fiberDim + foff+i] = stateVarsCell[varOffset+i];
      } // for
    } // for
  } // if/else

  topology::FieldBase::VectorFieldEnum multiType = topology::FieldBase::MULTI_OTHER;
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
      std::ostringstream msg;
      msg << "Bad vector field type '" << fieldType << "' for Material." << std::endl;
      throw std::logic_error(msg.str());
    } // switch
  field->vectorFieldType(multiType);

  PYLITH_METHOD_END;
} // getField
  
// ----------------------------------------------------------------------
// Get indices for physical property or state variable field.
void
pylith::materials::Material::_findField(int* propertyIndex,
					int* stateVarIndex,
					const char* name) const
{ // _findField
  PYLITH_METHOD_BEGIN;

  assert(propertyIndex);
  assert(stateVarIndex);

  *propertyIndex = -1;
  *stateVarIndex = -1;

  const std::string nameString = name;
  const int numProperties = _metadata.numProperties();
  for (int i=0; i < numProperties; ++i)
    if (nameString == _metadata.getProperty(i).name) {
      *propertyIndex = i;
      PYLITH_METHOD_END;
    } // if

  const int numStateVars = _metadata.numStateVars();
  for (int i=0; i < numStateVars; ++i)
    if (nameString == _metadata.getStateVar(i).name) {
      *stateVarIndex = i;
      PYLITH_METHOD_END;
    } // if

  PYLITH_METHOD_END;
} // _findField
  

// End of file 
