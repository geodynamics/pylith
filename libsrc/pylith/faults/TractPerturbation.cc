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
#include "pylith/topology/Fields.hh" // HOLDSA Fields
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
const pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >*
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
  _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);

  // Create section to hold time dependent values
  _parameters->add("value", "traction", topology::FieldBase::VERTICES_FIELD, spaceDim);
  _parameters->get("value").vectorFieldType(topology::FieldBase::VECTOR);
  _parameters->get("value").scale(pressureScale);
  _parameters->get("value").allocate();
  if (_dbInitial) {
    _parameters->add("initial", "traction_initial", topology::FieldBase::VERTICES_FIELD, spaceDim);
    _parameters->get("value").vectorFieldType(topology::FieldBase::VECTOR);
    _parameters->get("value").scale(pressureScale);
    _parameters->get("value").allocate();
  }
  if (_dbRate) {
    _parameters->add("rate", "traction_rate", topology::FieldBase::VERTICES_FIELD, spaceDim);
    _parameters->get("rate").vectorFieldType(topology::FieldBase::VECTOR);
    _parameters->get("rate").scale(rateScale);
    _parameters->get("rate").allocate();
    _parameters->add("rate time", "rate_start_time", topology::FieldBase::VERTICES_FIELD, 1);
    _parameters->get("rate time").vectorFieldType(topology::FieldBase::SCALAR);
    _parameters->get("rate time").scale(timeScale);
    _parameters->get("rate time").allocate();
  } // if
  if (_dbChange) {
    _parameters->add("change", "traction_change", topology::FieldBase::VERTICES_FIELD, spaceDim);
    _parameters->get("change").vectorFieldType(topology::FieldBase::VECTOR);
    _parameters->get("change").scale(pressureScale);
    _parameters->get("change").allocate();
    _parameters->add("change time", "change_start_time", topology::FieldBase::FACES_FIELD, 1);
    _parameters->get("change time").vectorFieldType(topology::FieldBase::SCALAR);
    _parameters->get("change time").scale(timeScale);
    _parameters->get("change time").allocate();
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
  DM             dmMesh = _parameters->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  PetscSection   initialSection, rateSection, rateTimeSection, changeSection, changeTimeSection;
  Vec            initialVec, rateVec, rateTimeVec, changeVec, changeTimeVec;
  PetscScalar   *valuesArray, *initialArray, *rateArray, *rateTimeArray, *changeArray, *changeTimeArray;

  PetscSection valuesSection     = _parameters->get("value").petscSection();
  Vec          valuesVec         = _parameters->get("value").localVector();
  assert(valuesSection);assert(valuesVec);
  err = VecGetArray(valuesVec, &valuesArray);CHECK_PETSC_ERROR(err);
  if (_dbInitial) {
    initialSection    = _parameters->get("initial").petscSection();
    initialVec        = _parameters->get("initial").localVector();
    assert(initialSection);assert(initialVec);
    err = VecGetArray(initialVec,    &initialArray);CHECK_PETSC_ERROR(err);
  }
  if (_dbRate) {
    rateSection       = _parameters->get("rate").petscSection();
    rateTimeSection   = _parameters->get("rate time").petscSection();
    rateVec           = _parameters->get("rate").localVector();
    rateTimeVec       = _parameters->get("rate time").localVector();
    assert(rateSection);assert(rateTimeSection);assert(rateVec);assert(rateTimeVec);
    err = VecGetArray(rateVec,       &rateArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(rateTimeVec,   &rateTimeArray);CHECK_PETSC_ERROR(err);
  }
  if (_dbChange) {
    changeSection     = _parameters->get("change").petscSection();
    changeTimeSection = _parameters->get("change time").petscSection();
    changeVec         = _parameters->get("change").localVector();
    changeTimeVec     = _parameters->get("change time").localVector();
    assert(changeSection);assert(changeTimeSection);assert(changeVec);assert(changeTimeVec);
    err = VecGetArray(changeVec,     &changeArray);CHECK_PETSC_ERROR(err);
    err = VecGetArray(changeTimeVec, &changeTimeArray);CHECK_PETSC_ERROR(err);
  }

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt vdof, voff;

    err = PetscSectionGetDof(valuesSection, v, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valuesSection, v, &voff);CHECK_PETSC_ERROR(err);
    assert(vdof == spaceDim);
    for(PetscInt d = 0; d < vdof; ++d)
      valuesArray[voff+d] = 0.0;

    // Contribution from initial value
    if (_dbInitial) {
      PetscInt idof, ioff;

      err = PetscSectionGetDof(initialSection, v, &idof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(initialSection, v, &ioff);CHECK_PETSC_ERROR(err);
      assert(idof == vdof);
      for(PetscInt d = 0; d < idof; ++d)
        valuesArray[voff+d] += initialArray[ioff+d];
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      PetscInt rdof, roff, rtdof, rtoff;

      err = PetscSectionGetDof(rateSection, v, &rdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateSection, v, &roff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(rateTimeSection, v, &rtdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateTimeSection, v, &rtoff);CHECK_PETSC_ERROR(err);
      assert(rdof == vdof);
      assert(rtdof == 1);

      const PylithScalar tRel = t - rateTimeArray[rtoff];
      if (tRel > 0.0)  // rate of change integrated over time
        for(PetscInt d = 0; d < spaceDim; ++d)
          valuesArray[voff+d] += rateArray[roff+d] * tRel;
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      PetscInt cdof, coff, ctdof, ctoff;

      err = PetscSectionGetDof(changeSection, v, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeSection, v, &coff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(changeTimeSection, v, &ctdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeTimeSection, v, &ctoff);CHECK_PETSC_ERROR(err);
      assert(cdof == vdof);
      assert(ctdof == 1);

      const PylithScalar tRel = t - changeTimeArray[ctoff];
      if (tRel >= 0) { // change in value over time
        PylithScalar scale = 1.0;
        if (0 != _dbTimeHistory) {
          PylithScalar tDim = tRel * timeScale;
          /* _getNormalizer().dimensionalize(&tDim, 1, timeScale);*/
          const int err = _dbTimeHistory->query(&scale, tDim);
          if (0 != err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < spaceDim; ++d)
          valuesArray[voff+d] += changeArray[coff+d] * scale;
      } // if
    } // if
  } // for
  err = VecRestoreArray(valuesVec,     &valuesArray);CHECK_PETSC_ERROR(err);
  if (_dbInitial) {
    err = VecRestoreArray(initialVec,    &initialArray);CHECK_PETSC_ERROR(err);
  }
  if (_dbRate) {
    err = VecRestoreArray(rateVec,       &rateArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(rateTimeVec,   &rateTimeArray);CHECK_PETSC_ERROR(err);
  }
  if (_dbChange) {
    err = VecRestoreArray(changeVec,     &changeArray);CHECK_PETSC_ERROR(err);
    err = VecRestoreArray(changeTimeVec, &changeTimeArray);CHECK_PETSC_ERROR(err);
  }
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
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::TractPerturbation::vertexField(const char* name,
					       const topology::SolutionFields* const fields)
{ // vertexField
  assert(_parameters);
  assert(name);

  if (0 == strcasecmp(name, "traction_initial_value"))
    return _parameters->get("initial");

  else if (0 == strcasecmp(name, "traction_rate_of_change"))
    return _parameters->get("rate");

  else if (0 == strcasecmp(name, "traction_change_in_value"))
    return _parameters->get("change");

  else if (0 == strcasecmp(name, "traction_rate_start_time"))
    return _parameters->get("rate time");

  else if (0 == strcasecmp(name, "traction_change_start_time"))
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
  DM             dmMesh = _parameters->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  const spatialdata::geocoords::CoordSys* cs = _parameters->mesh().coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();

  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  scalar_array valuesVertex(querySize);

  // Get sections.
  scalar_array coordsVertexGlobal(spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coordArray;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  PetscSection valueSection = _parameters->get("value").petscSection();
  Vec          valueVec     = _parameters->get("value").localVector();
  PetscScalar *tractionsArray;
  assert(valueSection);assert(valueVec);

  // Loop over cells in boundary mesh and perform queries.
  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(valueVec, &tractionsArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt cdof, coff;

    err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == cdof);
    for(PetscInt d = 0; d < cdof; ++d) {
      coordsVertexGlobal[d] = coordArray[coff+d];
    }
    normalizer.dimensionalize(&coordsVertexGlobal[0], coordsVertexGlobal.size(), lengthScale);
    
    valuesVertex = 0.0;
    err = db->query(&valuesVertex[0], querySize, &coordsVertexGlobal[0], spaceDim, cs);
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
    PetscInt vdof, voff;

    err = PetscSectionGetDof(valueSection, v, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valueSection, v, &voff);CHECK_PETSC_ERROR(err);
    assert(vdof == spaceDim);
    for(PetscInt d = 0; d < vdof; ++d)
      tractionsArray[voff+d] = valuesVertex[d];
  } // for
  err = VecRestoreArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
} // _queryDB

// End of file
