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

#include "TractPerturbation.hh" // implementation of object methods

#include "FaultCohesiveLagrange.hh" // USES faultToGlobal()

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

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
  PYLITH_METHOD_BEGIN;

  TimeDependent::deallocate();

  delete _parameters; _parameters = 0;

  PYLITH_METHOD_END;
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
const pylith::topology::Fields*
pylith::faults::TractPerturbation::parameterFields(void) const
{ // parameterFields
  return _parameters;
} // parameterFields

// ----------------------------------------------------------------------
// Initialize traction perturbation function.
void
pylith::faults::TractPerturbation::initialize(const topology::Mesh& faultMesh,
					      const topology::Field& faultOrientation, 
					      const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  PYLITH_METHOD_BEGIN;

  const PylithScalar pressureScale = normalizer.pressureScale();
  const PylithScalar timeScale = normalizer.timeScale();
  const PylithScalar rateScale = pressureScale / timeScale;
  _timeScale = timeScale;

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  delete _parameters; _parameters = new topology::Fields(faultMesh);assert(_parameters);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", topology::FieldBase::VERTICES_FIELD, spaceDim);
  topology::Field& value = _parameters->get("value");
  value.vectorFieldType(topology::FieldBase::VECTOR);
  value.scale(pressureScale);
  value.allocate();
  if (_dbInitial) {
    _parameters->add("initial", "traction_initial", topology::FieldBase::VERTICES_FIELD, spaceDim);
    topology::Field& initial = _parameters->get("initial");
    initial.vectorFieldType(topology::FieldBase::VECTOR);
    initial.scale(pressureScale);
    initial.allocate();
  }
  if (_dbRate) {
    _parameters->add("rate", "traction_rate", topology::FieldBase::VERTICES_FIELD, spaceDim);
    topology::Field& rate = _parameters->get("rate");
    rate.vectorFieldType(topology::FieldBase::VECTOR);
    rate.scale(rateScale);
    rate.allocate();
    _parameters->add("rate time", "rate_start_time", topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field& rateTime = _parameters->get("rate time");
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);
    rateTime.scale(timeScale);
    rateTime.allocate();
  } // if
  if (_dbChange) {
    _parameters->add("change", "traction_change", topology::FieldBase::VERTICES_FIELD, spaceDim);
    topology::Field& change = _parameters->get("change");
    change.vectorFieldType(topology::FieldBase::VECTOR);
    change.scale(pressureScale);
    change.allocate();
    _parameters->add("change time", "change_start_time", topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field& changeTime = _parameters->get("change time");
    changeTime.vectorFieldType(topology::FieldBase::SCALAR);
    changeTime.scale(timeScale);
    changeTime.allocate();
  } // if

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
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << " in TractPerturbation'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("initial", _dbInitial, spaceDim, pressureScale, normalizer);
    _dbInitial->close();
    pylith::topology::Field& initial = _parameters->get("initial");
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
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << " in TractPerturbation'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("rate", _dbRate, spaceDim, rateScale, normalizer);
    
    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale, normalizer);
    _dbRate->close();
    topology::Field& rate = _parameters->get("rate");
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
        std::ostringstream msg;
        msg << "Bad spatial dimension '" << spaceDim << " in TractPerturbation'." << std::endl;
        throw std::logic_error(msg.str());
      } // switch
    _queryDB("change", _dbChange, spaceDim, pressureScale, normalizer);
    
    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale, normalizer);
    _dbChange->close();
    topology::Field& change = _parameters->get("change");
    FaultCohesiveLagrange::faultToGlobal(&change, faultOrientation);

    if (_dbTimeHistory)
      _dbTimeHistory->open();
  } // if

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of Submesh.
void
pylith::faults::TractPerturbation::calculate(const PylithScalar t)
{ // calculate
  PYLITH_METHOD_BEGIN;

  assert(_parameters);

  // Get vertices.
  PetscDM dmMesh = _parameters->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();
  
  topology::Field* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::VecVisitorMesh* initialVisitor = (initialField) ? new topology::VecVisitorMesh(*initialField) : 0;
  PetscScalar* initialArray = (initialVisitor) ? initialVisitor->localArray() : NULL;

  topology::Field* rateField = (_dbRate) ? &_parameters->get("rate") : 0;
  topology::VecVisitorMesh* rateVisitor = (rateField) ? new topology::VecVisitorMesh(*rateField) : 0;
  PetscScalar* rateArray = (rateVisitor) ? rateVisitor->localArray() : NULL;

  topology::Field* rateTimeField = (_dbRate) ? &_parameters->get("rate time") : 0;
  topology::VecVisitorMesh* rateTimeVisitor = (rateTimeField) ? new topology::VecVisitorMesh(*rateTimeField) : 0;
  PetscScalar* rateTimeArray = (rateTimeVisitor) ? rateTimeVisitor->localArray() : NULL;

  topology::Field* changeField = (_dbChange) ? &_parameters->get("change") : 0;
  topology::VecVisitorMesh* changeVisitor = (changeField) ? new topology::VecVisitorMesh(*changeField) : 0;
  PetscScalar* changeArray = (changeVisitor) ? changeVisitor->localArray() : NULL;

  topology::Field* changeTimeField = (_dbChange) ? &_parameters->get("change time") : 0;
  topology::VecVisitorMesh* changeTimeVisitor = (changeTimeField) ? new topology::VecVisitorMesh(*changeTimeField) : 0;
  PetscScalar* changeTimeArray = (changeTimeVisitor) ? changeTimeVisitor->localArray() : NULL;

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt voff = valueVisitor.sectionOffset(v);
    assert(spaceDim == valueVisitor.sectionDof(v));
    for (PetscInt d = 0; d < spaceDim; ++d) {
      valueArray[voff+d] = 0.0;
    } // for

    // Contribution from initial value
    if (_dbInitial) {
      assert(initialVisitor);
      const PetscInt ioff = initialVisitor->sectionOffset(v);
      assert(spaceDim == initialVisitor->sectionDof(v));
      for (PetscInt d = 0; d < spaceDim; ++d) {
        valueArray[voff+d] += initialArray[ioff+d];
      } // for
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateVisitor);
      const PetscInt roff = rateVisitor->sectionOffset(v);
      assert(spaceDim == rateVisitor->sectionDof(v));
      assert(rateTimeVisitor);
      const PetscInt rtoff = rateTimeVisitor->sectionOffset(v);
      assert(1 == rateTimeVisitor->sectionDof(v));

      const PylithScalar tRel = t - rateTimeArray[rtoff];
      if (tRel > 0.0)  // rate of change integrated over time
	for(int iDim = 0; iDim < spaceDim; ++iDim) {
	  valueArray[voff+iDim] += rateArray[roff+iDim] * tRel;
	} // for
    } // if

    // Contribution from change of value
    if (_dbChange) {
      assert(changeVisitor);
      const PetscInt coff = changeVisitor->sectionOffset(v);
      assert(spaceDim == changeVisitor->sectionDof(v));
      const PetscInt ctoff = changeTimeVisitor->sectionOffset(v);
      assert(1 == changeTimeVisitor->sectionDof(v));

      const PylithScalar tRel = t - changeTimeArray[ctoff];
      if (tRel >= 0) { // change in value over time
	PylithScalar scale = 1.0;
	if (_dbTimeHistory) {
	  PylithScalar tDim = tRel*_timeScale;
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (err) {
	    std::ostringstream msg;
	    msg << "Error querying for time '" << tDim 
		<< "' in time history database '"
		<< _dbTimeHistory->label() << "'.";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	for (int iDim = 0; iDim < spaceDim; ++iDim) {
	  valueArray[voff+iDim] += changeArray[coff+iDim]*scale;
	} // for
      } // if
    } // if
  } // for

  delete initialVisitor; initialVisitor = 0;
  delete rateVisitor; rateVisitor = 0;
  delete rateTimeVisitor; rateTimeVisitor = 0;
  delete changeVisitor; changeVisitor = 0;
  delete changeTimeVisitor; changeTimeVisitor = 0;

  PYLITH_METHOD_END;
}  // calculate


// ----------------------------------------------------------------------
// Determine if perturbation has a given parameter.
bool
pylith::faults::TractPerturbation::hasParameter(const char* name) const
{ // hasParameter
  if (0 == strcasecmp(name, "traction_initial_value"))
    return (0 != _dbInitial);
  else if (0 == strcasecmp(name, "traction_rate_of_change"))
    return (0 != _dbRate);
  else if (0 == strcasecmp(name, "traction_change_in_value"))
    return (0 != _dbChange);
  else if (0 == strcasecmp(name, "traction_rate_start_time"))
    return (0 != _dbRate);
  else if (0 == strcasecmp(name, "traction_change_start_time"))
    return (0 != _dbChange);
  else
    return false;
} // hasParameter

// ----------------------------------------------------------------------
// Get vertex field with traction perturbation information.
const pylith::topology::Field&
pylith::faults::TractPerturbation::vertexField(const char* name,
					       const topology::SolutionFields* const fields)
{ // vertexField
  PYLITH_METHOD_BEGIN;

  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "traction_initial_value"))
    PYLITH_METHOD_RETURN(_parameters->get("initial"));

  else if (0 == strcasecmp(name, "traction_rate_of_change"))
    PYLITH_METHOD_RETURN(_parameters->get("rate"));

  else if (0 == strcasecmp(name, "traction_change_in_value"))
    PYLITH_METHOD_RETURN(_parameters->get("change"));

  else if (0 == strcasecmp(name, "traction_rate_start_time"))
    PYLITH_METHOD_RETURN(_parameters->get("rate time"));

  else if (0 == strcasecmp(name, "traction_change_start_time"))
    PYLITH_METHOD_RETURN(_parameters->get("change time"));

  else {
    std::ostringstream msg;
    msg << "Unknown field '" << name << "' requested for fault traction perturbation '" 
	<< _label << "'.";
    throw std::runtime_error(msg.str());
  } // else

  PYLITH_METHOD_RETURN(_parameters->get("traction")); // Satisfy method definition
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
  PYLITH_METHOD_BEGIN;

  assert(name);
  assert(db);
  assert(_parameters);

  // Get vertices.
  PetscDM dmMesh = _parameters->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();

  scalar_array coordsVertex(spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar *coordArray = coordsVisitor.localArray();

  topology::Field& parametersField = _parameters->get(name);
  topology::VecVisitorMesh parametersVisitor(parametersField);
  PetscScalar* parametersArray = parametersVisitor.localArray();

  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valueVertex(querySize);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    // Get dimensionalized coordinates of vertex
    const int coff = coordsVisitor.sectionOffset(v);
    assert(spaceDim == coordsVisitor.sectionDof(v));
    for (PetscInt d = 0; d < spaceDim; ++d) {
      coordsVertex[d] = coordArray[coff+d];
    } // for
    normalizer.dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);

    valueVertex = 0.0;
    int err = db->query(&valueVertex[0], valueVertex.size(), &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Error querying for '" << name << "' at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << coordsVertex[i];
      msg << ") using spatial database '" << db->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&valueVertex[0], valueVertex.size(), scale);

    // Update section
    const PetscInt off = parametersVisitor.sectionOffset(v);
    assert(querySize == parametersVisitor.sectionDof(v));
    for(int i = 0; i < querySize; ++i) {
      parametersArray[off+i] = valueVertex[i];
    } // for
  } // for

  PYLITH_METHOD_END;
} // _queryDB

// End of file
