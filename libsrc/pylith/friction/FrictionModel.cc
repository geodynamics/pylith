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

#include "FrictionModel.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES scalar_array, std::vector
#include "pylith/faults/FaultCohesiveLagrange.hh" // USES isClampedVertex()

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Default constructor.
pylith::friction::FrictionModel::FrictionModel(const materials::Metadata& metadata) :
  _dt(0.0),
  _normalizer(new spatialdata::units::Nondimensional),
  _metadata(metadata),
  _label(""),
  _dbProperties(0),
  _dbInitialState(0),
  _fieldsPropsStateVars(0),
  _propsFiberDim(0),
  _varsFiberDim(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::friction::FrictionModel::~FrictionModel(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::friction::FrictionModel::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _normalizer; _normalizer = 0;
  delete _fieldsPropsStateVars; _fieldsPropsStateVars = 0;
  _propsFiberDim = 0;
  _varsFiberDim = 0;

  _dbProperties = 0; // :TODO: Use shared pointer.
  _dbInitialState = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::friction::FrictionModel::normalizer(const spatialdata::units::Nondimensional& dim)
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
pylith::friction::FrictionModel::initialize(const topology::Mesh& faultMesh,
					    feassemble::Quadrature* quadrature)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_dbProperties);

  // Get vertices associated with friction interface
  PetscDM faultDMMesh = faultMesh.dmMesh();assert(faultDMMesh);
  topology::Stratum verticesStratum(faultDMMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  scalar_array coordsVertexGlobal(spaceDim);
  topology::CoordsVisitor coordsVisitor(faultDMMesh);
  PetscScalar* coordArray = coordsVisitor.localArray();

  // Query database for properties

  // Create fields to hold physical properties and state variables.
  delete _fieldsPropsStateVars; _fieldsPropsStateVars = new topology::Fields(faultMesh);assert(_fieldsPropsStateVars);
  _setupPropsStateVars();

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  scalar_array propertiesDBQuery(numDBProperties);
  scalar_array propertiesVertex(_propsFiberDim);

  // Setup database for querying for physical properties
  assert(_dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt coff = coordsVisitor.sectionOffset(v);
    assert(spaceDim == coordsVisitor.sectionDof(v));
    for (PetscInt d = 0; d < spaceDim; ++d) {
      coordsVertexGlobal[d] = coordArray[coff+d];
    } // for
    _normalizer->dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);

    int err = _dbProperties->query(&propertiesDBQuery[0], 
				   propertiesDBQuery.size(),
				   &coordsVertexGlobal[0], spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find parameters for physical properties at " << "(";
      for (int i = 0; i < spaceDim; ++i)
        msg << "  " << coordsVertexGlobal[i];
      msg << ") in friction model '" << _label << "' using spatial database '" << _dbProperties->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    assert(propertiesVertex.size() == propertiesDBQuery.size());
    _dbToProperties(&propertiesVertex[0], propertiesDBQuery);

    _nondimProperties(&propertiesVertex[0], propertiesVertex.size());
    PetscInt iOff = 0;

    for (int i=0; i < _metadata.numProperties(); ++i) {
      const materials::Metadata::ParamDescription& property = _metadata.getProperty(i);
      // TODO This needs to be an integer instead of a string
      topology::Field& propertyField = _fieldsPropsStateVars->get(property.name.c_str());
      topology::VecVisitorMesh propertyVisitor(propertyField);
      PetscScalar* propertyArray = propertyVisitor.localArray();
      const PetscInt off = propertyVisitor.sectionOffset(v);
      const PetscInt dof = propertyVisitor.sectionDof(v);
      for(PetscInt d = 0; d < dof; ++d, ++iOff) {
        propertyArray[off+d] += propertiesVertex[iOff];
      } // for
    } // for
  } // for
  // Close properties database
  _dbProperties->close();

  // Query database for initial state variables
  if (_dbInitialState) {

    // Create arrays for querying    
    const int numDBStateVars = _metadata.numDBStateVars();assert(numDBStateVars > 0);
    assert(_varsFiberDim > 0);
    scalar_array stateVarsDBQuery(numDBStateVars);
    scalar_array stateVarsVertex(_varsFiberDim);
    
    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(), _metadata.numDBStateVars());
    
    PetscDMLabel clamped = NULL;
    PetscErrorCode err = DMPlexGetLabel(faultDMMesh, "clamped", &clamped);PYLITH_CHECK_ERROR(err);

    for(PetscInt v = vStart; v < vEnd; ++v) {
      if (faults::FaultCohesiveLagrange::isClampedVertex(clamped, v)) {
	continue;
      } // if

      const PetscInt coff = coordsVisitor.sectionOffset(v);
      assert(spaceDim == coordsVisitor.sectionDof(v));
      for (PetscInt d = 0; d < spaceDim; ++d) {
	coordsVertexGlobal[d] = coordArray[coff+d];
      } // for
      _normalizer->dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);
      
      int err = _dbInitialState->query(&stateVarsDBQuery[0], numDBStateVars, &coordsVertexGlobal[0], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find initial state variables at " << "(";
        for (int i = 0; i < spaceDim; ++i)
          msg << "  " << coordsVertexGlobal[i];
        msg << ") in friction model '" << _label << "' using spatial database '" << _dbInitialState->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
      _dbToStateVars(&stateVarsVertex[0], stateVarsDBQuery);
      _nondimStateVars(&stateVarsVertex[0], stateVarsVertex.size());
      PetscInt iOff = 0;

      for (int i=0; i < _metadata.numStateVars(); ++i) {
        const materials::Metadata::ParamDescription& stateVar = _metadata.getStateVar(i);
        // TODO This needs to be an integer instead of a string
        topology::Field& stateVarField = _fieldsPropsStateVars->get(stateVar.name.c_str());
	topology::VecVisitorMesh stateVarVisitor(stateVarField);
	PetscScalar* stateVarArray = stateVarVisitor.localArray();
	const PetscInt off = stateVarVisitor.sectionOffset(v);
	const PetscInt dof = stateVarVisitor.sectionDof(v);
        for(PetscInt d = 0; d < dof; ++d, ++iOff) {
          stateVarArray[off+d] += stateVarsVertex[iOff];
        } // for
      } // for
    } // for
    // Close database
    _dbInitialState->close();
  } else if (_metadata.numDBStateVars()) {
    std::cerr << "WARNING: No initial state given for friction model '" << label() << "'. Using default value of zero." << std::endl;
  } // if/else

  // Setup buffers for restrict/update of properties and state variables.
  _propsStateVarsVertex.resize(_propsFiberDim+_varsFiberDim);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get the field with all properties and state variables.
const pylith::topology::Fields&
pylith::friction::FrictionModel::fieldsPropsStateVars(void) const
{ // fieldsPropsStateVars
  PYLITH_METHOD_BEGIN;

  assert(_fieldsPropsStateVars);
  PYLITH_METHOD_RETURN(*_fieldsPropsStateVars);
} // fieldsPropsStateVars

// ----------------------------------------------------------------------
// Check whether material has a field as a property.
bool
pylith::friction::FrictionModel::hasPropStateVar(const char* name)
{ // hasPropStateVar
  PYLITH_METHOD_BEGIN;

  if (_fieldsPropsStateVars) {
    return _fieldsPropsStateVars->hasField(name);
  } else {
    const std::string nameString = name;
    const int numProperties = _metadata.numProperties();
    for (int i=0; i < numProperties; ++i)
      if (_metadata.getProperty(i).name == nameString)
	PYLITH_METHOD_RETURN(true);
    const int numStateVars = _metadata.numStateVars();
    for (int i=0; i < numStateVars; ++i)
      if (_metadata.getStateVar(i).name == nameString)
	PYLITH_METHOD_RETURN(true);
  } // if/else

  PYLITH_METHOD_RETURN(false);
} // hasPropStateVar

// ----------------------------------------------------------------------
// Get metadta for physical properties or state variables.
const pylith::materials::Metadata&
pylith::friction::FrictionModel::getMetadata()
{ // getMetadata
  return _metadata;
} // getMetadata
  
// ----------------------------------------------------------------------
// Get physical property or state variable field.
const pylith::topology::Field&
pylith::friction::FrictionModel::getField(const char* name)
{ // getField
  PYLITH_METHOD_BEGIN;

  assert(name);
  assert(_fieldsPropsStateVars);

  PYLITH_METHOD_RETURN(_fieldsPropsStateVars->get(name));
} // getField
  
// ----------------------------------------------------------------------
// Retrieve properties and state variables for a point.
void
pylith::friction::FrictionModel::retrievePropsStateVars(const int point)
{ // retrievePropsStateVars
  PYLITH_METHOD_BEGIN;

  assert(_fieldsPropsStateVars);
  PetscInt iOff = 0;

  for (int i=0; i < _metadata.numProperties(); ++i) {
    const materials::Metadata::ParamDescription& property = _metadata.getProperty(i);
    // TODO This needs to be an integer instead of a string
    topology::Field& propertyField = _fieldsPropsStateVars->get(property.name.c_str());
    topology::VecVisitorMesh propertyVisitor(propertyField);
    PetscScalar* propertyArray = propertyVisitor.localArray();
    const PetscInt off = propertyVisitor.sectionOffset(point);
    const PetscInt dof = propertyVisitor.sectionDof(point);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      _propsStateVarsVertex[iOff] = propertyArray[off+d];
    } // for
  } // for
  for (int i=0; i < _metadata.numStateVars(); ++i) {
    const materials::Metadata::ParamDescription& stateVar = _metadata.getStateVar(i);
    // TODO This needs to be an integer instead of a string
    topology::Field& stateVarField = _fieldsPropsStateVars->get(stateVar.name.c_str());
    topology::VecVisitorMesh stateVarVisitor(stateVarField);
    PetscScalar* stateVarArray = stateVarVisitor.localArray();
    const PetscInt off = stateVarVisitor.sectionOffset(point);
    const PetscInt dof = stateVarVisitor.sectionDof(point);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      _propsStateVarsVertex[iOff] = stateVarArray[off+d];
    } // for
  } // for
  assert(_propsStateVarsVertex.size() == size_t(iOff));

  PYLITH_METHOD_END;
} // retrievePropsStateVars

// ----------------------------------------------------------------------
// Compute friction at vertex.
PylithScalar
pylith::friction::FrictionModel::calcFriction(const PylithScalar t,
					      const PylithScalar slip,
                                              const PylithScalar slipRate,
                                              const PylithScalar normalTraction)
{ // calcFriction
  PYLITH_METHOD_BEGIN;

  assert(_fieldsPropsStateVars);

  assert(size_t(_propsFiberDim+_varsFiberDim) == _propsStateVarsVertex.size());
  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  const PylithScalar* stateVarsVertex = (_varsFiberDim > 0) ?
    &_propsStateVarsVertex[_propsFiberDim] : 0;

  const PylithScalar friction = _calcFriction(t, slip, slipRate, normalTraction,
					      propertiesVertex, _propsFiberDim,
					      stateVarsVertex, _varsFiberDim);
  
  PYLITH_METHOD_RETURN(friction);
} // calcFriction

// ----------------------------------------------------------------------
// Compute derivative of friction with slip at vertex.
PylithScalar
pylith::friction::FrictionModel::calcFrictionDeriv(const PylithScalar t,
						   const PylithScalar slip,
						   const PylithScalar slipRate,
						   const PylithScalar normalTraction)
{ // calcFrictionDeriv
  PYLITH_METHOD_BEGIN;

  assert(_fieldsPropsStateVars);

  assert(size_t(_propsFiberDim+_varsFiberDim) == _propsStateVarsVertex.size());
  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  const PylithScalar* stateVarsVertex = (_varsFiberDim > 0) ?
    &_propsStateVarsVertex[_propsFiberDim] : 0;

  const PylithScalar frictionDeriv = _calcFrictionDeriv(t, slip, slipRate, normalTraction,
							propertiesVertex, _propsFiberDim,
							stateVarsVertex, _varsFiberDim);
  
  PYLITH_METHOD_RETURN(frictionDeriv);
} // calcFrictionDeriv

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::updateStateVars(const PylithScalar t,
						 const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const int vertex)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  assert(_fieldsPropsStateVars);
  if (0 == _varsFiberDim)
    PYLITH_METHOD_END;

  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  PylithScalar* stateVarsVertex = &_propsStateVarsVertex[_propsFiberDim];
  
  _updateStateVars(t, slip, slipRate, normalTraction,
		   &stateVarsVertex[0], _varsFiberDim,
		   &propertiesVertex[0], _propsFiberDim);

  PetscInt iOff = 0;

  for (int i=0; i < _metadata.numProperties(); ++i) {
    const materials::Metadata::ParamDescription& property = _metadata.getProperty(i);
    // TODO This needs to be an integer instead of a string
    topology::Field& propertyField = _fieldsPropsStateVars->get(property.name.c_str());
    topology::VecVisitorMesh propertyVisitor(propertyField);
    PetscScalar* propertyArray = propertyVisitor.localArray();
    const PetscInt off = propertyVisitor.sectionOffset(vertex);
    const PetscInt dof = propertyVisitor.sectionDof(vertex);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      propertyArray[off+d] = _propsStateVarsVertex[iOff];
    } // for
  } // for
  for (int i=0; i < _metadata.numStateVars(); ++i) {
    const materials::Metadata::ParamDescription& stateVar = _metadata.getStateVar(i);
    // TODO This needs to be an integer instead of a string
    topology::Field& stateVarField = _fieldsPropsStateVars->get(stateVar.name.c_str());
    topology::VecVisitorMesh stateVarVisitor(stateVarField);
    PetscScalar* stateVarArray = stateVarVisitor.localArray();
    const PetscInt off = stateVarVisitor.sectionOffset(vertex);
    const PetscInt dof = stateVarVisitor.sectionDof(vertex);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      stateVarArray[off+d] = _propsStateVarsVertex[iOff];
    } // for
  } // for
  assert(_propsStateVarsVertex.size() == size_t(iOff));

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::_updateStateVars(
				const PylithScalar t,
				const PylithScalar slip,
				const PylithScalar slipRate,
				const PylithScalar normalTraction,
				PylithScalar* const stateVars,
				const int numStateVars,
				const PylithScalar* properties,
				const int numProperties)
{ // _updateStateVars
} // _updateStateVars

// ----------------------------------------------------------------------
// Setup fields for physical properties and state variables.
void
pylith::friction::FrictionModel::_setupPropsStateVars(void)
{ // _setupPropsStateVars
  PYLITH_METHOD_BEGIN;

  // Determine number of values needed to store physical properties.
  const int numProperties = _metadata.numProperties();
  _propsFiberDim = 0;
  for (int i=0; i < numProperties; ++i)
    _propsFiberDim += _metadata.getProperty(i).fiberDim;
  assert(_propsFiberDim >= 0);
  
  // Determine scales for each physical property.
  scalar_array propertiesVertex(_propsFiberDim);
  for (int i=0; i < _propsFiberDim; ++i)
    propertiesVertex[i] = 1.0;
  _dimProperties(&propertiesVertex[0], propertiesVertex.size());

  // Determine number of values needed to store state variables.
  const int numStateVars = _metadata.numStateVars();
  _varsFiberDim = 0;
  for (int i=0; i < numStateVars; ++i)
    _varsFiberDim += _metadata.getStateVar(i).fiberDim;
  assert(_varsFiberDim >= 0);
  
  // Determine scales for each state variable.
  scalar_array stateVarsVertex(_varsFiberDim);
  for (int i=0; i < _varsFiberDim; ++i)
    stateVarsVertex[i] = 1.0;
  _dimStateVars(&stateVarsVertex[0], stateVarsVertex.size());

  // Setup fields
  assert(_fieldsPropsStateVars);

  for (int i=0, iScale=0; i < numProperties; ++i) {
    const materials::Metadata::ParamDescription& property = 
      _metadata.getProperty(i);
    _fieldsPropsStateVars->add(property.name.c_str(), property.name.c_str());
    topology::Field& propertyField = _fieldsPropsStateVars->get(property.name.c_str());
    propertyField.newSection(topology::FieldBase::VERTICES_FIELD, property.fiberDim);
    propertyField.allocate();
    propertyField.vectorFieldType(property.fieldType);
    propertyField.scale(propertiesVertex[iScale]);
    propertyField.zeroAll();
    iScale += property.fiberDim;
  } // for
  
  for (int i=0, iScale=0; i < numStateVars; ++i) {
    const materials::Metadata::ParamDescription& stateVar = 
      _metadata.getStateVar(i);
    _fieldsPropsStateVars->add(stateVar.name.c_str(), stateVar.name.c_str());
    topology::Field& stateVarField = _fieldsPropsStateVars->get(stateVar.name.c_str());
    stateVarField.newSection(topology::FieldBase::VERTICES_FIELD, stateVar.fiberDim);
    stateVarField.allocate();
    stateVarField.vectorFieldType(stateVar.fieldType);
    stateVarField.scale(stateVarsVertex[iScale]);
    stateVarField.zeroAll();
    iScale += stateVar.fiberDim;
  } // for
  assert(_varsFiberDim >= 0);

  PYLITH_METHOD_END;
} // _setupPropsStateVars


// End of file 
