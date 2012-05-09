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

#include "TractPerturbation.hh" // implementation of object methods

#include "FaultCohesiveLagrange.hh" // USES faultToGlobal()

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/FieldsNew.hh" // HOLDSA FieldsNew
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealUniformSection SubRealUniformSection;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::TractPerturbation::TractPerturbation(void) :
  _parameters(0),
  _timeScale(1.0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::TractPerturbation::~TractPerturbation(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::TractPerturbation::deallocate(void)
{ // deallocate
  TimeDependent::deallocate();

  delete _parameters; _parameters = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Set label for traction perturbation.
void
pylith::faults::TractPerturbation::label(const char* value)
{ // label
  _label = value;
} // label

// ----------------------------------------------------------------------
// Get parameter fields.
const pylith::topology::FieldsNew<pylith::topology::SubMesh>*
pylith::faults::TractPerturbation::parameterFields(void) const
{ // parameterFields
  return _parameters;
} // parameterFields

// ----------------------------------------------------------------------
// Initialize traction perturbation function.
void
pylith::faults::TractPerturbation::initialize(const topology::SubMesh& faultMesh,
					      const topology::Field<topology::SubMesh>& faultOrientation, 
					      const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  const PylithScalar pressureScale = normalizer.pressureScale();
  const PylithScalar timeScale = normalizer.timeScale();
  const PylithScalar rateScale = pressureScale / timeScale;
  _timeScale = timeScale;

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  delete _parameters; 
  _parameters = new topology::FieldsNew<topology::SubMesh>(faultMesh);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", spaceDim, topology::FieldBase::VECTOR, pressureScale);
  if (_dbInitial) 
    _parameters->add("initial", "traction_initial", spaceDim, topology::FieldBase::VECTOR, pressureScale);
  if (_dbRate) {
    _parameters->add("rate", "traction_rate", spaceDim, topology::FieldBase::VECTOR, rateScale);
    _parameters->add("rate time", "rate_start_time", 1, topology::FieldBase::SCALAR, timeScale);
  } // if
  if (_dbChange) {
    _parameters->add("change", "traction_change", spaceDim, topology::FieldBase::VECTOR, pressureScale);
    _parameters->add("change time", "change_start_time", 1, topology::FieldBase::SCALAR, timeScale);
  } // if
  _parameters->allocate(topology::FieldBase::VERTICES_FIELD, 0);
  const ALE::Obj<SubRealUniformSection>& parametersSection = _parameters->section();
  assert(!parametersSection.isNull());

  if (_dbInitial) { // Setup initial values, if provided.
    _dbInitial->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbInitial->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbInitial->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-leftlateral",
				    "traction-shear-updip",
				    "traction-normal"};
	_dbInitial->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in TractPerturbation.");
      } // switch
    _queryDB("initial", _dbInitial, spaceDim, pressureScale, normalizer);
    _dbInitial->close();
    pylith::topology::Field<pylith::topology::SubMesh>& initial = _parameters->get("initial");
    FaultCohesiveLagrange::faultToGlobal(&initial, faultOrientation);
  } // if

  if (_dbRate) { // Setup rate of change of values, if provided.
    _dbRate->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-rate-normal"};
	_dbRate->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-rate-shear", 
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-rate-shear-leftlateral",
				    "traction-rate-shear-updip",
				    "traction-rate-normal"};
	_dbRate->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in TractPerturbation.");
      } // switch
    _queryDB("rate", _dbRate, spaceDim, rateScale, normalizer);
    
    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale, normalizer);
    _dbRate->close();
    pylith::topology::Field<pylith::topology::SubMesh>& rate = _parameters->get("rate");
    FaultCohesiveLagrange::faultToGlobal(&rate, faultOrientation);
  } // if

  if (_dbChange) { // Setup change of values, if provided.
    _dbChange->open();
    switch (spaceDim)
      { // switch
      case 1 : {
	const char* valueNames[] = {"traction-normal"};
	_dbChange->queryVals(valueNames, 1);
	break;
      } // case 1
      case 2 : {
	const char* valueNames[] = {"traction-shear", "traction-normal"};
	_dbChange->queryVals(valueNames, 2);
	break;
      } // case 2
      case 3 : {
	const char* valueNames[] = {"traction-shear-leftlateral",
				    "traction-shear-updip",
				    "traction-normal"};
	_dbChange->queryVals(valueNames, 3);
	break;
      } // case 3
      default :
	std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad spatial dimension in TractPerturbation.");
      } // switch
    _queryDB("change", _dbChange, spaceDim, pressureScale, normalizer);
    
    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale, normalizer);
    _dbChange->close();
    pylith::topology::Field<pylith::topology::SubMesh>& change = _parameters->get("change");
    FaultCohesiveLagrange::faultToGlobal(&change, faultOrientation);

    if (_dbTimeHistory)
      _dbTimeHistory->open();
  } // if
  
  logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::faults::TractPerturbation::calculate(const PylithScalar t)
{ // calculate
  assert(_parameters);

  const PylithScalar timeScale = _timeScale;

  // Get vertices.
  const ALE::Obj<SieveSubMesh>& sieveMesh = _parameters->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const ALE::Obj<SubRealUniformSection>& parametersSection = _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  scalar_array parametersVertex(parametersFiberDim);
  
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(spaceDim == valueFiberDim);

  const int initialIndex = (_dbInitial) ? _parameters->sectionIndex("initial") : -1;
  const int initialFiberDim = (_dbInitial) ? _parameters->sectionFiberDim("initial") : 0;

  const int rateIndex = (_dbRate) ? _parameters->sectionIndex("rate") : -1;
  const int rateFiberDim = (_dbRate) ? _parameters->sectionFiberDim("rate") : 0;
  const int rateTimeIndex = (_dbRate) ? _parameters->sectionIndex("rate time") : -1;
  const int rateTimeFiberDim = (_dbRate) ? _parameters->sectionFiberDim("rate time") : 0;

  const int changeIndex = (_dbChange) ? _parameters->sectionIndex("change") : -1;
  const int changeFiberDim = (_dbChange) ? _parameters->sectionFiberDim("change") : 0;
  const int changeTimeIndex = (_dbChange) ? _parameters->sectionIndex("change time") : -1;
  const int changeTimeFiberDim = (_dbChange) ? _parameters->sectionFiberDim("change time") : 0;

  for(SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter != verticesEnd; ++v_iter) {
    assert(parametersFiberDim == parametersSection->getFiberDimension(*v_iter));
    parametersSection->restrictPoint(*v_iter, &parametersVertex[0], parametersVertex.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersVertex[valueIndex+i] = 0.0;

    // Contribution from initial value
    if (_dbInitial) {
      assert(initialIndex >= 0);
      assert(initialFiberDim == valueFiberDim);
      for (int i=0; i < initialFiberDim; ++i)
	parametersVertex[valueIndex+i] += parametersVertex[initialIndex+i];
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateIndex >= 0);
      assert(rateFiberDim == valueFiberDim);
      assert(rateTimeIndex >= 0);
      assert(rateTimeFiberDim == 1);
      
      const PylithScalar tRel = t - parametersVertex[rateTimeIndex];
      if (tRel > 0.0)  // rate of change integrated over time
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  parametersVertex[valueIndex+iDim] += parametersVertex[rateIndex+iDim] * tRel;
	} // for
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      assert(changeIndex >= 0);
      assert(changeFiberDim == valueFiberDim);
      assert(changeTimeIndex >= 0);
      assert(changeTimeFiberDim == 1);

      const PylithScalar tRel = t - parametersVertex[changeTimeIndex];
      if (tRel >= 0) { // change in value over time
	PylithScalar scale = 1.0;
	if (_dbTimeHistory) {
	  PylithScalar tDim = tRel*timeScale;
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (err) {
	    std::ostringstream msg;
	    msg << "Error querying for time '" << tDim 
		<< "' in time history database "
		<< _dbTimeHistory->label() << ".";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  parametersVertex[valueIndex+iDim] += parametersVertex[changeIndex+iDim] * scale;
	} // for
      } // if
    } // if
    
    parametersSection->updatePoint(*v_iter, &parametersVertex[0]);
  } // for
}  // calculate


// ----------------------------------------------------------------------
// Determine if perturbation has a given parameter.
bool
pylith::faults::TractPerturbation::hasParameter(const char* name) const
{ // hasParameter
  if (0 == strcasecmp(name, "initial_value"))
    return (0 != _dbInitial);
  else if (0 == strcasecmp(name, "rate_of_change"))
    return (0 != _dbRate);
  else if (0 == strcasecmp(name, "change_in_value"))
    return (0 != _dbChange);
  else if (0 == strcasecmp(name, "rate_start_time"))
    return (0 != _dbRate);
  else if (0 == strcasecmp(name, "change_start_time"))
    return (0 != _dbChange);
  else
    return false;
} // hasParameter

// ----------------------------------------------------------------------
// Get vertex field with traction perturbation information.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::TractPerturbation::vertexField(const char* name,
					       const topology::SolutionFields* const fields)
{ // vertexField
  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "initial_value"))
    return _parameters->get("initial");

  else if (0 == strcasecmp(name, "rate_of_change"))
    return _parameters->get("rate");

  else if (0 == strcasecmp(name, "change_in_value"))
    return _parameters->get("change");

  else if (0 == strcasecmp(name, "rate_start_time"))
    return _parameters->get("rate time");

  else if (0 == strcasecmp(name, "change_start_time"))
    return _parameters->get("change time");

  else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for fault traction perturbation '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  return _parameters->get("traction"); // Satisfy method definition
} // vertexField

// ----------------------------------------------------------------------
// Get label of boundary condition surface.
const char*
pylith::faults::TractPerturbation::_getLabel(void) const
{ // _getLabel
  return _label.c_str();
} // _getLabel

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::faults::TractPerturbation::_queryDB(const char* name,
					    spatialdata::spatialdb::SpatialDB* const db,
					    const int querySize,
					    const PylithScalar scale,
					    const spatialdata::units::Nondimensional& normalizer)
{ // _queryDB
  assert(name);
  assert(db);
  assert(_parameters);

  // Get vertices.
  const ALE::Obj<SieveSubMesh>& sieveMesh = _parameters->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveSubMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveSubMesh::label_sequence::iterator verticesEnd = vertices->end();

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();

  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valuesVertex(querySize);
  scalar_array coordsVertexGlobal(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& coordsSection = sieveMesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());

  const ALE::Obj<SubRealUniformSection>& parametersSection = _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex(name);
  const int valueFiberDim = _parameters->sectionFiberDim(name);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);
  scalar_array parametersVertex(parametersFiberDim);

  // Loop over cells in boundary mesh and perform queries.
  for (SieveSubMesh::label_sequence::iterator v_iter = verticesBegin; v_iter != verticesEnd; ++v_iter) {
    assert(spaceDim == coordsSection->getFiberDimension(*v_iter));
    coordsSection->restrictPoint(*v_iter, &coordsVertexGlobal[0], coordsVertexGlobal.size());
    normalizer.dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);
    
    valuesVertex = 0.0;
    const int err = db->query(&valuesVertex[0], querySize, &coordsVertexGlobal[0], spaceDim, cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find values at (";
      for (int i=0; i < spaceDim; ++i)
	msg << " " << coordsVertexGlobal[i];
      msg << ") for traction boundary condition " << _label << "\n"
	  << "using spatial database " << db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if

    normalizer.nondimensionalize(&valuesVertex[0], valuesVertex.size(), scale);

    // Update section
    assert(parametersFiberDim == parametersSection->getFiberDimension(*v_iter));
    parametersSection->restrictPoint(*v_iter, &parametersVertex[0], parametersVertex.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersVertex[valueIndex+i] = valuesVertex[i];
    
    parametersSection->updatePoint(*v_iter, &parametersVertex[0]);
  } // for
} // _queryDB

// End of file
