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
// Copyright (c) 2010 University of California, Davis
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
#include "pylith/utils/array.hh" // USES double_array, std::vector

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
			feassemble::Quadrature<topology::SubMesh>* quadrature,
			const topology::Field<topology::SubMesh>& area)
{ // initialize
  assert(0 != _dbProperties);
  assert(0 != quadrature);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  //logger.stagePush("Friction");

  // Get quadrature information
  const int numQuadPts = quadrature->numQuadPts();
  const int numBasis = quadrature->numBasis();
  const int spaceDim = quadrature->spaceDim();
  const double_array& quadWts = quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cells associated with friction interface
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = faultMesh.sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices =
    faultSieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin =
    vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates =
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
        coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // Create fields to hold physical properties and state variables.
  delete _fieldsPropsStateVars; 
  _fieldsPropsStateVars = new topology::FieldsNew<topology::SubMesh>(faultMesh);
  assert(_fieldsPropsStateVars);
  _setupPropsStateVars();
  assert(_fieldsPropsStateVars->fiberDim() > 0);
  _fieldsPropsStateVars->allocate(vertices);
  
  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();

  // Query database for properties

  // Create arrays for querying.
  const int numDBProperties = _metadata.numDBProperties();
  double_array propertiesDBQuery(numDBProperties);
  double_array propertiesQuadPt(_propsFiberDim);

  const int fieldsFiberDim = _fieldsPropsStateVars->fiberDim();
  assert(_propsFiberDim + _varsFiberDim == fieldsFiberDim);
  const ALE::Obj<SubRealUniformSection>& fieldsSection = 
    _fieldsPropsStateVars->section();
  assert(!fieldsSection.isNull());
  double_array fieldsCell(numBasis*fieldsFiberDim);
  UpdateAddVisitor fieldsVisitor(*fieldsSection, &fieldsCell[0]);

  // Setup database for querying for physical properties
  assert(_dbProperties);
  _dbProperties->open();
  _dbProperties->queryVals(_metadata.dbProperties(),
			   _metadata.numDBProperties());

  fieldsSection->zero();
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    fieldsCell = 0.0;

    const double_array& quadPtsNonDim = quadrature->quadPts();
    quadPtsGlobal = quadPtsNonDim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
        lengthScale);

    // Loop over quadrature points in cell and query database
    for (int iQuadPt=0, index=0;
        iQuadPt < numQuadPts;
        ++iQuadPt, index+=spaceDim) {
      int err = _dbProperties->query(&propertiesDBQuery[0], numDBProperties,
          &quadPtsGlobal[index], spaceDim, cs);
      if (err) {
        std::ostringstream msg;
        msg << "Could not find parameters for physical properties at \n" << "(";
        for (int i = 0; i < spaceDim; ++i)
          msg << "  " << quadPtsGlobal[index + i];
        msg << ") in friction model " << _label << "\n"
            << "using spatial database '" << _dbProperties->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if
      _dbToProperties(&propertiesQuadPt[0], propertiesDBQuery);

      // Get cell geometry information that depends on cell
      const double_array& basis = quadrature->basis();
      const double_array& jacobianDet = quadrature->jacobianDet();

      // Compute properties weighted by area
      const double wt = quadWts[iQuadPt] * jacobianDet[iQuadPt];
      for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
        const double dArea = wt * basis[iQuadPt*numBasis+iBasis];
        for (int iProp = 0; iProp < _propsFiberDim; ++iProp)
          fieldsCell[iBasis*fieldsFiberDim+iProp]
              += propertiesQuadPt[iProp] * dArea;
      } // for
    } // for
    fieldsVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, fieldsVisitor);
  } // for
  // Close properties database
  _dbProperties->close();

  // Query database for initial state variables
  if (_dbInitialState) {
    assert(_varsFiberDim > 0);

    // Create arrays for querying
    const int numDBStateVars = _metadata.numDBStateVars();
    double_array stateVarsDBQuery(numDBStateVars);
    double_array stateVarsQuadPt(_varsFiberDim);

    // Setup database for querying for initial state variables
    _dbInitialState->open();
    _dbInitialState->queryVals(_metadata.dbStateVars(),
			       _metadata.numDBStateVars());

    for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
        c_iter != cellsEnd;
        ++c_iter) {
      // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
      quadrature->retrieveGeometry(*c_iter);
#else
      coordsVisitor.clear();
      faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
      quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

      const double_array& quadPtsNonDim = quadrature->quadPts();
      quadPtsGlobal = quadPtsNonDim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      fieldsCell = 0.0;

      // Loop over quadrature points in cell and query database
      for (int iQuadPt=0, index=0;
          iQuadPt < numQuadPts;
          ++iQuadPt, index+=spaceDim) {
        int err = _dbInitialState->query(&stateVarsDBQuery[0], numDBStateVars,
					 &quadPtsGlobal[index], spaceDim, cs);
        if (err) {
          std::ostringstream msg;
          msg << "Could not find initial state variables at \n" << "(";
          for (int i = 0; i < spaceDim; ++i)
            msg << "  " << quadPtsGlobal[index + i];
          msg << ") in friction model " << _label << "\n"
              << "using spatial database '" << _dbInitialState->label() << "'.";
          throw std::runtime_error(msg.str());
        } // if
        _dbToStateVars(&stateVarsQuadPt[0], stateVarsDBQuery);

        // Get cell geometry information that depends on cell
        const double_array& basis = quadrature->basis();
        const double_array& jacobianDet = quadrature->jacobianDet();

        // Compute state variables weighted by area
        const double wt = quadWts[iQuadPt] * jacobianDet[iQuadPt];
        for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
          const double dArea = wt * basis[iQuadPt*numBasis+iBasis];
          for (int iVar = 0; iVar < numDBStateVars; ++iVar)
            fieldsCell[iBasis*fieldsFiberDim+_propsFiberDim+iVar]
                += stateVarsDBQuery[iVar] * dArea;
        } // for
      } // for
      fieldsVisitor.clear();
      faultSieveMesh->updateClosure(*c_iter, fieldsVisitor);
    } // for
    // Close database
    _dbInitialState->close();
  } // if

  _fieldsPropsStateVars->complete(); // Assemble contributions

  // Loop over vertices and divide by area to get weighted values,
  // nondimensionalize properties and state variables.
  const ALE::Obj<RealSection>& areaSection = area.section();
  assert(!areaSection.isNull());

  double_array fieldsVertex(fieldsFiberDim);
  double areaVertex = 0.0;
  for (SieveSubMesh::label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    fieldsSection->restrictPoint(*v_iter, 
				 &fieldsVertex[0], fieldsVertex.size());

    _nondimProperties(&fieldsVertex[0], _propsFiberDim);
    if (_varsFiberDim > 0)
      _nondimStateVars(&fieldsVertex[_propsFiberDim], _varsFiberDim);

    assert(1 == areaSection->getFiberDimension(*v_iter));
    const double* areaVertex = areaSection->restrictPoint(*v_iter);
    assert(*areaVertex > 0.0);
    fieldsVertex /= *areaVertex;

    assert(fieldsFiberDim == fieldsSection->getFiberDimension(*v_iter));
    fieldsSection->updatePoint(*v_iter, &fieldsVertex[0]);
  } // for

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
double
pylith::friction::FrictionModel::calcFriction(const double slip,
                                              const double slipRate,
                                              const double normalTraction)
{ // calcFriction
  assert(_fieldsPropsStateVars);

  assert(_propsFiberDim+_varsFiberDim == _propsStateVarsVertex.size());
  const double* propertiesVertex = &_propsStateVarsVertex[0];
  const double* stateVarsVertex = (_varsFiberDim > 0) ?
    &_propsStateVarsVertex[_propsFiberDim] : 0;

  const double friction =
    _calcFriction(slip, slipRate, normalTraction,
		  propertiesVertex, _propsFiberDim,
		  stateVarsVertex, _varsFiberDim);
  
  return friction;
} // calcFriction

// ----------------------------------------------------------------------
// Update state variables (for next time step).
void
pylith::friction::FrictionModel::updateStateVars(const double slip,
						 const double slipRate,
						 const double normalTraction,
						 const int vertex)
{ // updateStateVars
  assert(_fieldsPropsStateVars);

  if (0 == _varsFiberDim)
    return;

  const ALE::Obj<SubRealUniformSection>& fieldsSection =
    _fieldsPropsStateVars->section();
  assert(!fieldsSection.isNull());

  const double* propertiesVertex = &_propsStateVarsVertex[0];
  double* stateVarsVertex = &_propsStateVarsVertex[_propsFiberDim];
  
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
pylith::friction::FrictionModel::_updateStateVars(const double slip,
    const double slipRate,
    const double normalTraction,
    double* const stateVars,
    const int numStateVars,
    const double* properties,
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
  double_array propertiesVertex(_propsFiberDim);
  for (int i=0; i < _propsFiberDim; ++i)
    propertiesVertex[i] = 1.0;
  _nondimProperties(&propertiesVertex[0], propertiesVertex.size());

  // Determine number of values needed to store state variables.
  const int numStateVars = _metadata.numStateVars();
  _varsFiberDim = 0;
  for (int i=0; i < numStateVars; ++i)
    _varsFiberDim += _metadata.getStateVar(i).fiberDim;
  assert(_varsFiberDim >= 0);
  
  // Determine scales for each state variable.
  double_array stateVarsVertex(_varsFiberDim);
  for (int i=0; i < _varsFiberDim; ++i)
    stateVarsVertex[i] = 1.0;
  _nondimStateVars(&stateVarsVertex[0], stateVarsVertex.size());

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
