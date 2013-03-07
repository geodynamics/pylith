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

#include "FrictionModel.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES Mesh
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
typedef pylith::topology::Mesh::SieveSubMesh SieveSubMesh;

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
  delete _normalizer; _normalizer = 0;
  delete _fieldsPropsStateVars; _fieldsPropsStateVars = 0;
  _propsFiberDim = 0;
  _varsFiberDim = 0;

  _dbProperties = 0; // :TODO: Use shared pointer.
  _dbInitialState = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Set scales used to nondimensionalize physical properties.
void
pylith::friction::FrictionModel::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::friction::FrictionModel::initialize(
			const topology::SubMesh& faultMesh,
			feassemble::Quadrature<topology::SubMesh>* quadrature)
{ // initialize
  assert(0 != _dbProperties);
  PetscErrorCode err;
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.stagePush("Friction");

  // Get vertices associated with friction interface
  DM       faultDMMesh = faultMesh.dmMesh();
  PetscInt vStart, vEnd;

  err = DMPlexGetDepthStratum(faultDMMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  scalar_array coordsVertexGlobal(spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar* coordArray;
  err = DMPlexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(faultDMMesh, &coordVec);CHECK_PETSC_ERROR(err);

  // Query database for properties

  // Create fields to hold physical properties and state variables.
  delete _fieldsPropsStateVars; 
  _fieldsPropsStateVars = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(_fieldsPropsStateVars);
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

  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt coff;
    PetscInt cdof;

    err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == cdof);
    for (PetscInt d = 0; d < cdof; ++d) {
      coordsVertexGlobal[d] = coordArray[coff+d];
    } // for
    _normalizer->dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);

    int err = _dbProperties->query(&propertiesDBQuery[0], 
				   propertiesDBQuery.size(),
				   &coordsVertexGlobal[0], spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find parameters for physical properties at \n" << "(";
      for (int i = 0; i < spaceDim; ++i)
        msg << "  " << coordsVertexGlobal[i];
      msg << ") in friction model " << _label << "\n"
          << "using spatial database '" << _dbProperties->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    assert(propertiesVertex.size() == propertiesDBQuery.size());
    _dbToProperties(&propertiesVertex[0], propertiesDBQuery);

    _nondimProperties(&propertiesVertex[0], propertiesVertex.size());
    PetscInt iOff = 0;

    for (int i=0; i < _metadata.numProperties(); ++i) {
      const materials::Metadata::ParamDescription& property = 
        _metadata.getProperty(i);
      // TODO This needs to be an integer instead of a string
      topology::Field<topology::SubMesh>& prop = _fieldsPropsStateVars->get(property.name.c_str());
      PetscSection section = prop.petscSection();
      Vec          vec     = prop.localVector();
      PetscScalar *a;
      PetscInt     dof, off;

      err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
      err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
      for(PetscInt d = 0; d < dof; ++d, ++iOff) {
        a[off+d] += propertiesVertex[iOff];
      } // for
    } // for
  } // for
  // Close properties database
  _dbProperties->close();

  // Query database for initial state variables
  if (_dbInitialState) {
    assert(_varsFiberDim > 0);

    // Create arrays for querying
    const int numDBStateVars = _metadata.numDBStateVars();
    scalar_array stateVarsDBQuery(numDBStateVars);
    scalar_array stateVarsVertex(_varsFiberDim);
    
    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(),
			       _metadata.numDBStateVars());
    
    for(PetscInt v = vStart; v < vEnd; ++v) {
      PetscInt coff;
      PetscInt cdof;

      err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
      assert(spaceDim == cdof);
      for (PetscInt d = 0; d < cdof; ++d) {
	coordsVertexGlobal[d] = coordArray[coff+d];
      } // for
      _normalizer->dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);
      
      int err = _dbInitialState->query(&stateVarsDBQuery[0], numDBStateVars,
				       &coordsVertexGlobal[0], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find initial state variables at \n" << "(";
        for (int i = 0; i < spaceDim; ++i)
          msg << "  " << coordsVertexGlobal[i];
        msg << ") in friction model " << _label << "\n"
            << "using spatial database '" << _dbInitialState->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
      _dbToStateVars(&stateVarsVertex[0], stateVarsDBQuery);
      _nondimStateVars(&stateVarsVertex[0], stateVarsVertex.size());
      PetscInt iOff = 0;

      for (int i=0; i < _metadata.numStateVars(); ++i) {
        const materials::Metadata::ParamDescription& stateVar = 
          _metadata.getStateVar(i);
        // TODO This needs to be an integer instead of a string
        topology::Field<topology::SubMesh>& sv = _fieldsPropsStateVars->get(stateVar.name.c_str());
        PetscSection section = sv.petscSection();
        Vec          vec     = sv.localVector();
        PetscScalar *a;
        PetscInt     dof, off;
        
        err = PetscSectionGetDof(section, v, &dof);CHECK_PETSC_ERROR(err);
        err = PetscSectionGetOffset(section, v, &off);CHECK_PETSC_ERROR(err);
        err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
        for(PetscInt d = 0; d < dof; ++d, ++iOff) {
          a[off+d] += stateVarsVertex[iOff];
        }
      } // for
    } // for
    // Close database
    _dbInitialState->close();
  } else if (_metadata.numDBStateVars()) {
    std::cerr << "WARNING: No initial state given for friction model '" << label() << "'. Using default value of zero." << std::endl;
  } // if/else

  err = VecRestoreArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);

  // Setup buffers for restrict/update of properties and state variables.
  _propsStateVarsVertex.resize(_propsFiberDim+_varsFiberDim);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Get the field with all properties and state variables.
const pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >&
pylith::friction::FrictionModel::fieldsPropsStateVars(void) const
{ // fieldsPropsStateVars
  assert(_fieldsPropsStateVars);
  return *_fieldsPropsStateVars;
} // fieldsPropsStateVars

// ----------------------------------------------------------------------
// Check whether material has a field as a property.
bool
pylith::friction::FrictionModel::hasPropStateVar(const char* name)
{ // hasPropStateVar
  if (_fieldsPropsStateVars) {
    return _fieldsPropsStateVars->hasField(name);
  } else {
    const std::string nameString = name;
    const int numProperties = _metadata.numProperties();
    for (int i=0; i < numProperties; ++i)
      if (_metadata.getProperty(i).name == nameString)
	return true;
    const int numStateVars = _metadata.numStateVars();
    for (int i=0; i < numStateVars; ++i)
      if (_metadata.getStateVar(i).name == nameString)
	return true;
  } // if/else

  return false;
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
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::friction::FrictionModel::getField(const char* name)
{ // getField
  assert(name);
  assert(_fieldsPropsStateVars);

  return _fieldsPropsStateVars->get(name);
} // getField
  
// ----------------------------------------------------------------------
// Retrieve properties and state variables for a point.
void
pylith::friction::FrictionModel::retrievePropsStateVars(const int point)
{ // retrievePropsStateVars
  assert(_fieldsPropsStateVars);
  PetscErrorCode err;
  PetscInt iOff = 0;

  for (int i=0; i < _metadata.numProperties(); ++i) {
    const materials::Metadata::ParamDescription& property = 
      _metadata.getProperty(i);
    // TODO This needs to be an integer instead of a string
    topology::Field<topology::SubMesh>& prop = _fieldsPropsStateVars->get(property.name.c_str());
    PetscSection section = prop.petscSection();
    Vec          vec     = prop.localVector();
    PetscScalar *a;
    PetscInt     dof, off;

    err = PetscSectionGetDof(section, point, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, point, &off);CHECK_PETSC_ERROR(err);
    err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      _propsStateVarsVertex[iOff] = a[off+d];
    }
  }
  for (int i=0; i < _metadata.numStateVars(); ++i) {
    const materials::Metadata::ParamDescription& stateVar = 
      _metadata.getStateVar(i);
    // TODO This needs to be an integer instead of a string
    topology::Field<topology::SubMesh>& sv = _fieldsPropsStateVars->get(stateVar.name.c_str());
    PetscSection section = sv.petscSection();
    Vec          vec     = sv.localVector();
    PetscScalar *a;
    PetscInt     dof, off;

    err = PetscSectionGetDof(section, point, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, point, &off);CHECK_PETSC_ERROR(err);
    err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      _propsStateVarsVertex[iOff] = a[off+d];
    }
  } // for
  assert(_propsStateVarsVertex.size() == iOff);
} // retrievePropsStateVars

// ----------------------------------------------------------------------
// Compute friction at vertex.
PylithScalar
pylith::friction::FrictionModel::calcFriction(const PylithScalar t,
					      const PylithScalar slip,
                                              const PylithScalar slipRate,
                                              const PylithScalar normalTraction)
{ // calcFriction
  assert(_fieldsPropsStateVars);

  assert(_propsFiberDim+_varsFiberDim == _propsStateVarsVertex.size());
  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  const PylithScalar* stateVarsVertex = (_varsFiberDim > 0) ?
    &_propsStateVarsVertex[_propsFiberDim] : 0;

  const PylithScalar friction =
    _calcFriction(t, slip, slipRate, normalTraction,
		  propertiesVertex, _propsFiberDim,
		  stateVarsVertex, _varsFiberDim);
  
  return friction;
} // calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::updateStateVars(const PylithScalar t,
						 const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const int vertex)
{ // updateStateVars
  assert(_fieldsPropsStateVars);
  PetscErrorCode err;
  if (0 == _varsFiberDim)
    return;

  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  PylithScalar* stateVarsVertex = &_propsStateVarsVertex[_propsFiberDim];
  
  _updateStateVars(t, slip, slipRate, normalTraction,
		   &stateVarsVertex[0], _varsFiberDim,
		   &propertiesVertex[0], _propsFiberDim);

  PetscInt iOff = 0;

  for (int i=0; i < _metadata.numProperties(); ++i) {
    const materials::Metadata::ParamDescription& property = 
      _metadata.getProperty(i);
    // TODO This needs to be an integer instead of a string
    topology::Field<topology::SubMesh>& prop = _fieldsPropsStateVars->get(property.name.c_str());
    PetscSection section = prop.petscSection();
    Vec          vec     = prop.localVector();
    PetscScalar *a;
    PetscInt     dof, off;

    err = PetscSectionGetDof(section, vertex, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, vertex, &off);CHECK_PETSC_ERROR(err);
    err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      a[off+d] = _propsStateVarsVertex[iOff];
    }
  }
  for (int i=0; i < _metadata.numStateVars(); ++i) {
    const materials::Metadata::ParamDescription& stateVar = 
      _metadata.getStateVar(i);
    // TODO This needs to be an integer instead of a string
    topology::Field<topology::SubMesh>& sv = _fieldsPropsStateVars->get(stateVar.name.c_str());
    PetscSection section = sv.petscSection();
    Vec          vec     = sv.localVector();
    PetscScalar *a;
    PetscInt     dof, off;

    err = PetscSectionGetDof(section, vertex, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(section, vertex, &off);CHECK_PETSC_ERROR(err);
    err = VecGetArray(vec, &a);CHECK_PETSC_ERROR(err);
    for(PetscInt d = 0; d < dof; ++d, ++iOff) {
      a[off+d] = _propsStateVarsVertex[iOff];
    }
  } // for
  assert(_propsStateVarsVertex.size() == iOff);
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
    topology::Field<topology::SubMesh>& prop = _fieldsPropsStateVars->get(property.name.c_str());
    prop.newSection(topology::FieldBase::VERTICES_FIELD, property.fiberDim);
    prop.allocate();
    prop.vectorFieldType(property.fieldType);
    prop.scale(propertiesVertex[iScale]);
    prop.zero();
    iScale += property.fiberDim;
  } // for
  
  for (int i=0, iScale=0; i < numStateVars; ++i) {
    const materials::Metadata::ParamDescription& stateVar = 
      _metadata.getStateVar(i);
    _fieldsPropsStateVars->add(stateVar.name.c_str(), stateVar.name.c_str());
    topology::Field<topology::SubMesh>& sv = _fieldsPropsStateVars->get(stateVar.name.c_str());
    sv.newSection(topology::FieldBase::VERTICES_FIELD, stateVar.fiberDim);
    sv.allocate();
    sv.vectorFieldType(stateVar.fieldType);
    sv.scale(stateVarsVertex[iScale]);
    sv.zero();
    iScale += stateVar.fiberDim;
  } // for
  assert(_varsFiberDim >= 0);
} // _setupPropsStateVars


// End of file 
