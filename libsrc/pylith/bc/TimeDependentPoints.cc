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

#include "TimeDependentPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

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
  BoundaryConditionPoints::deallocate();
  TimeDependent::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set indices of vertices with point forces.
void
pylith::bc::TimeDependentPoints::bcDOF(const int* flags,
				       const int size)
{ // bcDOF
  if (size > 0)
    assert(0 != flags);

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
  const PylithScalar timeScale = _getNormalizer().timeScale();
  const PylithScalar rateScale = valueScale / timeScale;
  PetscErrorCode err = 0;

  const int numPoints = _points.size();
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

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  delete _parameters;
  _parameters = new topology::Fields<topology::Field<topology::Mesh> >(mesh);

  _parameters->add("value", "value", topology::FieldBase::VERTICES_FIELD, numBCDOF);
  topology::Field<topology::Mesh>& value = _parameters->get("value");
  value.vectorFieldType(topology::FieldBase::OTHER);
  value.scale(valueScale);
  value.allocate();
  PetscVec valueVec = value.localVector();assert(valueVec);
  err = VecSet(valueVec, 0.0);CHECK_PETSC_ERROR(err);
  
  if (_dbInitial) {
    const std::string& fieldLabel = std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field<topology::Mesh>& initial = _parameters->get("initial");
    initial.vectorFieldType(topology::FieldBase::OTHER);
    initial.scale(valueScale);
    initial.allocate();
    PetscVec initialVec = initial.localVector();assert(initialVec);
    err = VecSet(initialVec, 0.0);CHECK_PETSC_ERROR(err);
  } // if
  if (_dbRate) {
    const std::string& fieldLabel = std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field<topology::Mesh>& rate = _parameters->get("rate");
    rate.vectorFieldType(topology::FieldBase::OTHER);
    rate.scale(rateScale);
    rate.allocate();
    PetscVec rateVec = rate.localVector();assert(rateVec);
    err = VecSet(rateVec, 0.0);CHECK_PETSC_ERROR(err);
    const std::string& timeLabel = std::string("rate_time_") + std::string(fieldName);    
    _parameters->add("rate time", timeLabel.c_str(), topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field<topology::Mesh>& rateTime = _parameters->get("rate time");
    rateTime.vectorFieldType(topology::FieldBase::SCALAR);
    rateTime.scale(timeScale);
    rateTime.allocate();
    PetscVec rateTimeVec = rateTime.localVector();assert(rateTimeVec);
    err = VecSet(rateTimeVec, 0.0);CHECK_PETSC_ERROR(err);
  } // if
  if (_dbChange) {
    const std::string& fieldLabel = std::string("change_") + std::string(fieldName);
    _parameters->add("change", fieldLabel.c_str(), topology::FieldBase::VERTICES_FIELD, numBCDOF);
    topology::Field<topology::Mesh>& change = _parameters->get("change");
    change.vectorFieldType(topology::FieldBase::OTHER);
    change.scale(valueScale);
    change.allocate();
    PetscVec changeVec = change.localVector();assert(changeVec);
    err = VecSet(changeVec, 0.0);CHECK_PETSC_ERROR(err);
    const std::string& timeLabel = std::string("change_time_") + std::string(fieldName);
    _parameters->add("change time", timeLabel.c_str(), topology::FieldBase::VERTICES_FIELD, 1);
    topology::Field<topology::Mesh>& changeTime = _parameters->get("change time");
    changeTime.vectorFieldType(topology::FieldBase::SCALAR);
    changeTime.scale(timeScale);
    changeTime.allocate();
    PetscVec changeTimeVec = _parameters->get("change time").localVector();assert(changeTimeVec);
    err = VecSet(changeTimeVec, 0.0);CHECK_PETSC_ERROR(err);
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
    
    if (0 != _dbTimeHistory)
      _dbTimeHistory->open();
  } // if
  
  // Dellocate memory
  for (int i=0; i < numBCDOF; ++i) {
    delete[] valueNames[i]; valueNames[i] = 0;
    delete[] rateNames[i]; rateNames[i] = 0;
  } // for
  delete[] valueNames; valueNames = 0;
  delete[] rateNames; rateNames = 0;

  logger.stagePop();
} // _queryDatabases

// ----------------------------------------------------------------------
// Query database for values.
void
pylith::bc::TimeDependentPoints::_queryDB(const char* name,
				 spatialdata::spatialdb::SpatialDB* const db,
				 const int querySize,
				 const PylithScalar scale)
{ // _queryDB
  assert(name);
  assert(db);
  assert(_parameters);

  const topology::Mesh& mesh = _parameters->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = _getNormalizer().lengthScale();

  scalar_array coordsVertex(spaceDim);
  DM           dmMesh = mesh.dmMesh();
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coordArray;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);

  PetscSection parametersSection = _parameters->get(name).petscSection();
  Vec          parametersVec     = _parameters->get(name).localVector();
  assert(parametersSection);assert(parametersVec);
#if 0
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex(name);
  const int valueFiberDim = _parameters->sectionFiberDim(name);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);
  scalar_array parametersVertex(parametersFiberDim);
#endif

  scalar_array valueVertex(querySize);

  const int numPoints = _points.size();
  PetscScalar   *array;
  err = VecGetArray(parametersVec, &array);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    PetscInt cdof, coff;

    // Get dimensionalized coordinates of vertex
    err = PetscSectionGetDof(coordSection, _points[iPoint], &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, _points[iPoint], &coff);CHECK_PETSC_ERROR(err);
    assert(cdof == spaceDim);
    for (PetscInt d = 0; d < spaceDim; ++d) {coordsVertex[d] = coordArray[coff+d];}
    _getNormalizer().dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);
    int err = db->query(&valueVertex[0], valueVertex.size(), &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Error querying for '" << name << "' at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    _getNormalizer().nondimensionalize(&valueVertex[0], valueVertex.size(), scale);

    // Update section
    PetscInt dof, off;

    err = PetscSectionGetDof(parametersSection, _points[iPoint], &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(parametersSection, _points[iPoint], &off);CHECK_PETSC_ERROR(err);
    //assert(parametersFiberDim == dof);
    for(int i = 0; i < dof; ++i)
      array[off+i] = valueVertex[i];
  } // for
  err = VecRestoreArray(parametersVec, &array);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
} // _queryDB

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValue(const PylithScalar t)
{ // _calculateValue
  assert(_parameters);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const PylithScalar timeScale = _getNormalizer().timeScale();
  PetscSection   initialSection, rateSection, rateTimeSection, changeSection, changeTimeSection;
  Vec            initialVec, rateVec, rateTimeVec, changeVec, changeTimeVec;
  PetscScalar   *valuesArray, *initialArray, *rateArray, *rateTimeArray, *changeArray, *changeTimeArray;
  PetscErrorCode err;

  PetscSection valuesSection = _parameters->get("value").petscSection();
  Vec          valuesVec     = _parameters->get("value").localVector();
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

  for(int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label
    PetscInt vdof, voff;

    err = PetscSectionGetDof(valuesSection, p_bc, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valuesSection, p_bc, &voff);CHECK_PETSC_ERROR(err);
    // Contribution from initial value
    if (_dbInitial) {
      PetscInt idof, ioff;

      err = PetscSectionGetDof(initialSection, p_bc, &idof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(initialSection, p_bc, &ioff);CHECK_PETSC_ERROR(err);
      assert(idof == vdof);
      for(PetscInt d = 0; d < vdof; ++d)
        valuesArray[voff+d] += initialArray[ioff+d];
    } // if
    
    // Contribution from rate of change of value
    if (_dbRate) {
      PetscInt rdof, roff, rtdof, rtoff;

      err = PetscSectionGetDof(rateSection, p_bc, &rdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateSection, p_bc, &roff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(rateTimeSection, p_bc, &rtdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateTimeSection, p_bc, &rtoff);CHECK_PETSC_ERROR(err);
      assert(rdof == vdof);
      assert(rtdof == 1);
      const PylithScalar tRel = t - rateTimeArray[rtoff];
      if (tRel > 0.0)  // rate of change integrated over time
        for(PetscInt d = 0; d < vdof; ++d)
          valuesArray[voff+d] += rateArray[roff+d] * tRel;
    } // if

    // Contribution from change of value
    if (_dbChange) {
      PetscInt cdof, coff, ctdof, ctoff;

      err = PetscSectionGetDof(changeSection, p_bc, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeSection, p_bc, &coff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(changeTimeSection, p_bc, &ctdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeTimeSection, p_bc, &ctoff);CHECK_PETSC_ERROR(err);
      assert(cdof == vdof);
      assert(ctdof == 1);
      const PylithScalar tRel = t - changeTimeArray[ctoff];
      if (tRel >= 0) { // change in value over time
        PylithScalar scale = 1.0;
        if (0 != _dbTimeHistory) {
          PylithScalar tDim = tRel;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          const int err = _dbTimeHistory->query(&scale, tDim);
          if (0 != err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < vdof; ++d)
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
}  // _calculateValue

// ----------------------------------------------------------------------
// Calculate increment in temporal and spatial variation of value over
// the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValueIncr(const PylithScalar t0,
						     const PylithScalar t1)
{ // _calculateValueIncr
  assert(_parameters);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const PylithScalar timeScale = _getNormalizer().timeScale();
  PetscSection   initialSection, rateSection, rateTimeSection, changeSection, changeTimeSection;
  Vec            initialVec, rateVec, rateTimeVec, changeVec, changeTimeVec;
  PetscScalar   *valuesArray, *initialArray, *rateArray, *rateTimeArray, *changeArray, *changeTimeArray;
  PetscErrorCode err;

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

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label
    PetscInt vdof, voff;

    err = PetscSectionGetDof(valuesSection, p_bc, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valuesSection, p_bc, &voff);CHECK_PETSC_ERROR(err);
    // No contribution from initial value
    
    // Contribution from rate of change of value
    if (_dbRate) {
      PetscInt rdof, roff, rtdof, rtoff;

      err = PetscSectionGetDof(rateSection, p_bc, &rdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateSection, p_bc, &roff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(rateTimeSection, p_bc, &rtdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(rateTimeSection, p_bc, &rtoff);CHECK_PETSC_ERROR(err);
      assert(rdof == vdof);
      assert(rtdof == 1);
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
        for(PetscInt d = 0; d < vdof; ++d)
          valuesArray[voff+d] += rateArray[roff+d] * tIncr;
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      PetscInt cdof, coff, ctdof, ctoff;

      err = PetscSectionGetDof(changeSection, p_bc, &cdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeSection, p_bc, &coff);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(changeTimeSection, p_bc, &ctdof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(changeTimeSection, p_bc, &ctoff);CHECK_PETSC_ERROR(err);
      assert(cdof == vdof);
      assert(ctdof == 1);
      const PylithScalar tChange = changeTimeArray[ctoff];
      if (t0 >= tChange) { // increment is after change starts
        PylithScalar scale0 = 1.0;
        PylithScalar scale1 = 1.0;
        if (0 != _dbTimeHistory) {
          PylithScalar tDim = t0 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          int err = _dbTimeHistory->query(&scale0, tDim);
          if (0 != err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
          tDim = t1 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          err = _dbTimeHistory->query(&scale1, tDim);
          if (0 != err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < vdof; ++d)
          valuesArray[voff+d] += changeArray[coff+d] * (scale1 - scale0);
      } else if (t1 >= tChange) { // increment spans when change starts
        PylithScalar scale1 = 1.0;
        if (0 != _dbTimeHistory) {
          PylithScalar tDim = t1 - tChange;
          _getNormalizer().dimensionalize(&tDim, 1, timeScale);
          int err = _dbTimeHistory->query(&scale1, tDim);
          if (0 != err) {
            std::ostringstream msg;
            msg << "Error querying for time '" << tDim 
                << "' in time history database "
                << _dbTimeHistory->label() << ".";
            throw std::runtime_error(msg.str());
          } // if
        } // if
        for(PetscInt d = 0; d < vdof; ++d)
          valuesArray[voff+d] += changeArray[coff+d] * scale1;
      } // if/else
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
}  // _calculateValueIncr


// End of file 
