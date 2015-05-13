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

#include "TimeDependentPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::TimeDependentPoints::TimeDependentPoints(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::TimeDependentPoints::~TimeDependentPoints(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::TimeDependentPoints::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  BoundaryConditionPoints::deallocate();
  TimeDependent::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set indices of vertices with point forces.
void
pylith::bc::TimeDependentPoints::bcDOF(const int* flags,
				       const int size)
{ // bcDOF
  if (size > 0)
    assert(flags);

  _bcDOF.resize(size);
  for (int i=0; i < size; ++i)
    _bcDOF[i] = flags[i];
} // bcDOF

// ----------------------------------------------------------------------
// Query databases for parameters.
void
pylith::bc::TimeDependentPoints::_queryDatabases(const topology::Mesh& mesh,
						 const PylithScalar valueScale,
						 const char* fieldName)
{ // _queryDatabases
  PYLITH_METHOD_BEGIN;

  const PylithScalar timeScale = _getNormalizer().timeScale();
  const PylithScalar rateScale = valueScale / timeScale;

  const int numBCDOF = _bcDOF.size();
  char** valueNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  char** rateNames = (numBCDOF > 0) ? new char*[numBCDOF] : 0;
  const std::string& valuePrefix = std::string(fieldName) + "-";
  const std::string& ratePrefix = std::string(fieldName) + "-rate-";
  for (int i=0; i < numBCDOF; ++i) {
    std::string suffix;
    switch (_bcDOF[i])
      { // switch
      case 0 :
	suffix = "x";
	break;
      case 1 :
	suffix = "y";
	break;
      case 2 :
	suffix = "z";
	break;
      } // switch
    std::string name = valuePrefix + suffix;
    int size = 1 + name.length();
    valueNames[i] = new char[size];
    strcpy(valueNames[i], name.c_str());

    name = ratePrefix + suffix;
    size = 1 + name.length();
    rateNames[i] = new char[size];
    strcpy(rateNames[i], name.c_str());
  } // for

  delete _parameters; _parameters = new topology::Fields(mesh);assert(_parameters);

  PetscErrorCode err = 0;

  _parameters->add("value", "value", topology::FieldBase::VERTICES_FIELD, numBCDOF);
  topology::Field& value = _parameters->get("value");
  value.vectorFieldType(topology::FieldBase::OTHER);
  value.scale(valueScale);
  value.allocate();
  PetscVec valueVec = value.localVector();assert(valueVec);
  err = VecSet(valueVec, 0.0);PYLITH_CHECK_ERROR(err);
  
  if (_dbInitial) {
    const std::string& fieldLabel = std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field& initial = _parameters->get("initial");
    initial.vectorFieldType(topology::FieldBase::OTHER);
    initial.scale(valueScale);
    initial.allocate();
    PetscVec initialVec = initial.localVector();assert(initialVec);
    err = VecSet(initialVec, 0.0);PYLITH_CHECK_ERROR(err);
  } // if
  if (_dbRate) {
    const std::string& fieldLabel = std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field& rate = _parameters->get("rate");
    rate.vectorFieldType(topology::FieldBase::OTHER);
    rate.scale(rateScale);
    rate.allocate();
    PetscVec rateVec = rate.localVector();assert(rateVec);
    err = VecSet(rateVec, 0.0);PYLITH_CHECK_ERROR(err);
    const std::string& timeLabel = std::string("rate_time_") + std::string(fieldName);    
    _parameters->add("rate time", timeLabel.c_str(), topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field& rateTime = _parameters->get("rate time");
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);
    rateTime.scale(timeScale);
    rateTime.allocate();
    PetscVec rateTimeVec = rateTime.localVector();assert(rateTimeVec);
    err = VecSet(rateTimeVec, 0.0);PYLITH_CHECK_ERROR(err);
  } // if
  if (_dbChange) {
    const std::string& fieldLabel = std::string("change_") + std::string(fieldName);
    _parameters->add("change", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field& change = _parameters->get("change");
    change.vectorFieldType(topology::FieldBase::OTHER);
    change.scale(valueScale);
    change.allocate();
    PetscVec changeVec = change.localVector();assert(changeVec);
    err = VecSet(changeVec, 0.0);PYLITH_CHECK_ERROR(err);
    const std::string& timeLabel = std::string("change_time_") + std::string(fieldName);
    _parameters->add("change time", timeLabel.c_str(), topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field& changeTime = _parameters->get("change time");
    changeTime.vectorFieldType(topology::FieldBase::SCALAR);
    changeTime.scale(timeScale);
    changeTime.allocate();
    PetscVec changeTimeVec = _parameters->get("change time").localVector();assert(changeTimeVec);
    err = VecSet(changeTimeVec, 0.0);PYLITH_CHECK_ERROR(err);
  } // if
  
  if (_dbInitial) { // Setup initial values, if provided.
    _dbInitial->open();
    _dbInitial->queryVals(valueNames, numBCDOF);
    _queryDB("initial", _dbInitial, numBCDOF, valueScale);
    _dbInitial->close();
  } // if

  if (_dbRate) { // Setup rate of change of values, if provided.
    _dbRate->open();
    _dbRate->queryVals(rateNames, numBCDOF);
    _queryDB("rate", _dbRate, numBCDOF, rateScale);
    
    const char* timeNames[1] = { "rate-start-time" };
    _dbRate->queryVals(timeNames, 1);
    _queryDB("rate time", _dbRate, 1, timeScale);
    _dbRate->close();
  } // if
  
  if (_dbChange) { // Setup change of values, if provided.
    _dbChange->open();
    _dbChange->queryVals(valueNames, numBCDOF);
    _queryDB("change", _dbChange, numBCDOF, valueScale);
    
    const char* timeNames[1] = { "change-start-time" };
    _dbChange->queryVals(timeNames, 1);
    _queryDB("change time", _dbChange, 1, timeScale);
    _dbChange->close();
    
    if (_dbTimeHistory)
      _dbTimeHistory->open();
  } // if
  
  // Dellocate memory
  for (int i=0; i < numBCDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
    delete[] rateNames[i]; rateNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;
  delete[] rateNames; rateNames = 0;

  PYLITH_METHOD_END;
} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::TimeDependentPoints::_queryDB(const char* name,
					  spatialdata::spatialdb::SpatialDB* const db,
					  const int querySize,
					  const PylithScalar scale)
{ // _queryDB
  PYLITH_METHOD_BEGIN;

  assert(name);
  assert(db);
  assert(_parameters);

  const topology::Mesh& mesh = _parameters->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  const spatialdata::units::Nondimensional& normalizer = _getNormalizer();
  const PylithScalar lengthScale = normalizer.lengthScale();

  scalar_array coordsVertex(spaceDim);
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar *coordArray = coordsVisitor.localArray();

  topology::Field& parametersField = _parameters->get(name);
  topology::VecVisitorMesh parametersVisitor(parametersField);
  PetscScalar* parametersArray = parametersVisitor.localArray();

  scalar_array valueVertex(querySize);
  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get dimensionalized coordinates of vertex
    const int coff = coordsVisitor.sectionOffset(_points[iPoint]);
    assert(spaceDim == coordsVisitor.sectionDof(_points[iPoint]));
    for (PetscInt d = 0; d < spaceDim; ++d) {
      coordsVertex[d] = coordArray[coff+d];
    } // for
    normalizer.dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);

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
    const PetscInt off = parametersVisitor.sectionOffset(_points[iPoint]);
    assert(querySize == parametersVisitor.sectionDof(_points[iPoint]));
    for(int i = 0; i < querySize; ++i) {
      parametersArray[off+i] = valueVertex[i];
    } // for
  } // for

  PYLITH_METHOD_END;
} // _queryDB

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValue(const PylithScalar t)
{ // _calculateValue
  PYLITH_METHOD_BEGIN;

  assert(_parameters);

  const PylithScalar timeScale = _getNormalizer().timeScale();

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

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  for(int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label

    const PetscInt voff = valueVisitor.sectionOffset(p_bc);
    assert(numBCDOF == valueVisitor.sectionDof(p_bc));
    for (PetscInt d = 0; d < numBCDOF; ++d) {
      valueArray[voff+d] = 0.0;
    } // for

    // Contribution from initial value
    if (_dbInitial) {
      assert(initialVisitor);
      const PetscInt ioff = initialVisitor->sectionOffset(p_bc);
      assert(numBCDOF == initialVisitor->sectionDof(p_bc));
      for (PetscInt d = 0; d < numBCDOF; ++d) {
        valueArray[voff+d] += initialArray[ioff+d];
      } // for
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateVisitor);
      const PetscInt roff = rateVisitor->sectionOffset(p_bc);
      assert(numBCDOF == rateVisitor->sectionDof(p_bc));
      assert(rateTimeVisitor);
      const PetscInt rtoff = rateTimeVisitor->sectionOffset(p_bc);
      assert(1 == rateTimeVisitor->sectionDof(p_bc));

      const PylithScalar tRel = t - rateTimeArray[rtoff];
      if (tRel > 0.0)  // rate of change integrated over time
	for(int iDim = 0; iDim < numBCDOF; ++iDim) {
	  valueArray[voff+iDim] += rateArray[roff+iDim] * tRel;
	} // for
    } // if

    // Contribution from change of value
    if (_dbChange) {
      assert(changeVisitor);
      const PetscInt coff = changeVisitor->sectionOffset(p_bc);
      assert(numBCDOF == changeVisitor->sectionDof(p_bc));
      const PetscInt ctoff = changeTimeVisitor->sectionOffset(p_bc);
      assert(1 == changeTimeVisitor->sectionDof(p_bc));

      const PylithScalar tRel = t - changeTimeArray[ctoff];
      if (tRel >= 0) { // change in value over time
	PylithScalar scale = 1.0;
	if (_dbTimeHistory) {
	  PylithScalar tDim = tRel;
	  _getNormalizer().dimensionalize(&tDim, 1, timeScale);
	  const int err = _dbTimeHistory->query(&scale, tDim);
	  if (err) {
	    std::ostringstream msg;
	    msg << "Error querying for time '" << tDim 
		<< "' in time history database '"
		<< _dbTimeHistory->label() << "'.";
	    throw std::runtime_error(msg.str());
	  } // if
	} // if
	for (int iDim = 0; iDim < numBCDOF; ++iDim) {
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
}  // _calculateValue

// ----------------------------------------------------------------------
// Calculate increment in temporal and spatial variation of value over
// the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValueIncr(const PylithScalar t0,
						     const PylithScalar t1)
{ // _calculateValueIncr
  PYLITH_METHOD_BEGIN;

  assert(_parameters);

  const PylithScalar timeScale = _getNormalizer().timeScale();

  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();
  
  topology::Field* initialField = (_dbInitial) ? &_parameters->get("initial") : 0;
  topology::VecVisitorMesh* initialVisitor = (initialField) ? new topology::VecVisitorMesh(*initialField) : 0;

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

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label

    const PetscInt voff = valueVisitor.sectionOffset(p_bc);
    assert(numBCDOF == valueVisitor.sectionDof(p_bc));
    for (PetscInt d = 0; d < numBCDOF; ++d) {
      valueArray[voff+d] = 0.0;
    } // for

    // No contribution from initial value
    
    // Contribution from rate of change of value
    if (_dbRate) {
      const PetscInt roff = rateVisitor->sectionOffset(p_bc);
      assert(numBCDOF == rateVisitor->sectionDof(p_bc));
      assert(rateTimeVisitor);
      const PetscInt rtoff = rateTimeVisitor->sectionOffset(p_bc);
      assert(1 == rateTimeVisitor->sectionDof(p_bc));

      // Account for when rate dependence begins.
      const PylithScalar tRate = rateTimeArray[rtoff];
      PylithScalar tIncr = 0.0;
      if (t0 > tRate) // rate dependence for t0 to t1
        tIncr = t1 - t0;
      else if (t1 > tRate) // rate dependence for tRef to t1
        tIncr = t1 - tRate;
      else
        tIncr = 0.0; // no rate dependence for t0 to t1
      
      if (tIncr > 0.0)  // rate of change integrated over time
        for(PetscInt d = 0; d < numBCDOF; ++d)
          valueArray[voff+d] += rateArray[roff+d] * tIncr;
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      assert(changeVisitor);
      const PetscInt coff = changeVisitor->sectionOffset(p_bc);
      assert(numBCDOF == changeVisitor->sectionDof(p_bc));
      const PetscInt ctoff = changeTimeVisitor->sectionOffset(p_bc);
      assert(1 == changeTimeVisitor->sectionDof(p_bc));

      const PylithScalar tChange = changeTimeArray[ctoff];
      if (t0 >= tChange) { // increment is after change starts
        PylithScalar scale0 = 1.0;
        PylithScalar scale1 = 1.0;
        if (_dbTimeHistory) {
          PylithScalar tDim = t0 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          int err = _dbTimeHistory->query(&scale0, tDim);
          if (err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database '"
                << _dbTimeHistory->label() << "'.";
            throw std::runtime_error(msg.str());
          } // if
          tDim = t1 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          err = _dbTimeHistory->query(&scale1, tDim);
          if (err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database '"
                << _dbTimeHistory->label() << "'.";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < numBCDOF; ++d)
          valueArray[voff+d] += changeArray[coff+d] * (scale1 - scale0);
      } else if (t1 >= tChange) { // increment spans when change starts
        PylithScalar scale1 = 1.0;
        if (_dbTimeHistory) {
          PylithScalar tDim = t1 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          int err = _dbTimeHistory->query(&scale1, tDim);
          if (err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < numBCDOF; ++d)
          valueArray[voff+d] += changeArray[coff+d] * scale1;
      } // if/else
    } // if
  } // for

  delete initialVisitor; initialVisitor = 0;
  delete rateVisitor; rateVisitor = 0;
  delete rateTimeVisitor; rateTimeVisitor = 0;
  delete changeVisitor; changeVisitor = 0;
  delete changeTimeVisitor; changeTimeVisitor = 0;

  PYLITH_METHOD_END;
}  // _calculateValueIncr


// End of file 
