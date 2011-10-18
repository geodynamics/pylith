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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FrictionModel.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
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
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::FieldsNew<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;

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

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.stagePush("Friction");

  // Get vertices associated with friction interface
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  assert(0 != _normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  scalar_array coordsVertex(spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Query database for properties

  // Create fields to hold physical properties and state variables.
  delete _fieldsPropsStateVars; 
  _fieldsPropsStateVars = new topology::FieldsNew<topology::SubMesh>(faultMesh);
  assert(_fieldsPropsStateVars);
  _setupPropsStateVars();

  const int fieldsFiberDim = _fieldsPropsStateVars->fiberDim();
  assert(fieldsFiberDim > 0);
  _fieldsPropsStateVars->allocate(vertices);
  
  assert(_propsFiberDim + _varsFiberDim == fieldsFiberDim);
  const ALE::Obj<SubRealUniformSection>& fieldsSection = 
    _fieldsPropsStateVars->section();
  assert(!fieldsSection.isNull());
  scalar_array fieldsVertex(fieldsFiberDim);

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  scalar_array propertiesDBQuery(numDBProperties);
  scalar_array propertiesVertex(_propsFiberDim);

  // Setup database for querying for physical properties
  assert(_dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  fieldsSection->zero();
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, &coordsVertex[0], coordsVertex.size());
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(),
				lengthScale);


    int err = _dbProperties->query(&propertiesDBQuery[0], 
				   propertiesDBQuery.size(),
				   &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find parameters for physical properties at \n" << "(";
      for (int i = 0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") in friction model " << _label << "\n"
	  << "using spatial database '" << _dbProperties->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    assert(propertiesVertex.size() == propertiesDBQuery.size());
    _dbToProperties(&propertiesVertex[0], propertiesDBQuery);

    _nondimProperties(&propertiesVertex[0], propertiesVertex.size());

    fieldsVertex = 0.0;
    for (int iProp=0; iProp < _propsFiberDim; ++iProp)
      fieldsVertex[iProp] = propertiesVertex[iProp];
    assert(fieldsVertex.size() == fieldsSection->getFiberDimension(*v_iter));
    fieldsSection->updateAddPoint(*v_iter, &fieldsVertex[0]);
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
    
    for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
        v_iter != verticesEnd;
        ++v_iter) {
      coordinates->restrictPoint(*v_iter, &coordsVertex[0], 
				 coordsVertex.size());
      _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(),
				  lengthScale);
      
      int err = _dbInitialState->query(&stateVarsDBQuery[0], numDBStateVars,
				       &coordsVertex[0], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find initial state variables at \n" << "(";
	for (int i = 0; i < spaceDim; ++i)
	  msg << "  " << coordsVertex[i];
	msg << ") in friction model " << _label << "\n"
	    << "using spatial database '" << _dbInitialState->label() << "'.";
	throw std::runtime_error(msg.str());
      } // if
      _dbToStateVars(&stateVarsVertex[0], stateVarsDBQuery);
      _nondimStateVars(&stateVarsVertex[0], stateVarsVertex.size());

      fieldsVertex = 0.0;
      for (int iVar=0; iVar < _varsFiberDim; ++iVar)
	fieldsVertex[_propsFiberDim+iVar] = stateVarsVertex[iVar];
      assert(fieldsVertex.size() == fieldsSection->getFiberDimension(*v_iter));
      fieldsSection->updateAddPoint(*v_iter, &fieldsVertex[0]);
    } // for
    // Close database
    _dbInitialState->close();
  } // if

  // Setup buffers for restrict/update of properties and state variables.
  _propsStateVarsVertex.resize(fieldsFiberDim);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Get the field with all properties and state variables.
const pylith::topology::FieldsNew<pylith::topology::SubMesh>&
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

  const ALE::Obj<SubRealUniformSection>& fieldsSection =
    _fieldsPropsStateVars->section();
  assert(!fieldsSection.isNull());

  assert(_propsStateVarsVertex.size() ==
	 fieldsSection->getFiberDimension(point));
  fieldsSection->restrictPoint(point, &_propsStateVarsVertex[0],
			       _propsStateVarsVertex.size());
} // retrievePropsStateVars

// ----------------------------------------------------------------------
// Compute friction at vertex.
PylithScalar
pylith::friction::FrictionModel::calcFriction(const PylithScalar slip,
                                              const PylithScalar slipRate,
                                              const PylithScalar normalTraction)
{ // calcFriction
  assert(_fieldsPropsStateVars);

  assert(_propsFiberDim+_varsFiberDim == _propsStateVarsVertex.size());
  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  const PylithScalar* stateVarsVertex = (_varsFiberDim > 0) ?
    &_propsStateVarsVertex[_propsFiberDim] : 0;

  const PylithScalar friction =
    _calcFriction(slip, slipRate, normalTraction,
		  propertiesVertex, _propsFiberDim,
		  stateVarsVertex, _varsFiberDim);
  
  return friction;
} // calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::updateStateVars(const PylithScalar slip,
						 const PylithScalar slipRate,
						 const PylithScalar normalTraction,
						 const int vertex)
{ // updateStateVars
  assert(_fieldsPropsStateVars);

  if (0 == _varsFiberDim)
    return;

  const ALE::Obj<SubRealUniformSection>& fieldsSection =
    _fieldsPropsStateVars->section();
  assert(!fieldsSection.isNull());

  const PylithScalar* propertiesVertex = &_propsStateVarsVertex[0];
  PylithScalar* stateVarsVertex = &_propsStateVarsVertex[_propsFiberDim];
  
  _updateStateVars(slip, slipRate, normalTraction,
		   &stateVarsVertex[0], _varsFiberDim,
		   &propertiesVertex[0], _propsFiberDim);

  assert(_propsStateVarsVertex.size() == 
	 fieldsSection->getFiberDimension(vertex));
  fieldsSection->updatePoint(vertex, &_propsStateVarsVertex[0]);
} // updateStateVars

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::_updateStateVars(const PylithScalar slip,
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
    _fieldsPropsStateVars->add(property.name.c_str(), property.name.c_str(),
			       property.fiberDim, property.fieldType,
			       propertiesVertex[iScale]);
    iScale += property.fiberDim;
  } // for
  
  for (int i=0, iScale=0; i < numStateVars; ++i) {
    const materials::Metadata::ParamDescription& stateVar = 
      _metadata.getStateVar(i);
    _fieldsPropsStateVars->add(stateVar.name.c_str(), stateVar.name.c_str(),
			       stateVar.fiberDim, stateVar.fieldType,
			       stateVarsVertex[iScale]);
    iScale += stateVar.fiberDim;
  } // for
  assert(_varsFiberDim >= 0);
} // _setupPropsStateVars


// End of file 
