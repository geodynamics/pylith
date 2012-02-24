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
#include "pylith/topology/FieldsNew.hh" // USES FieldsNew
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealUniformSection RealUniformSection;
typedef pylith::topology::Mesh::RealSection RealSection;

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
						 const double valueScale,
						 const char* fieldName)
{ // _queryDatabases
  const double timeScale = _getNormalizer().timeScale();
  const double rateScale = valueScale / timeScale;

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
  _parameters = new topology::FieldsNew<topology::Mesh>(mesh);

  _parameters->add("value", "value", numBCDOF, topology::FieldBase::OTHER,
		   valueScale);
  
  if (_dbInitial) {
    const std::string& fieldLabel =
      std::string("initial_") + std::string(fieldName);
    _parameters->add("initial", fieldLabel.c_str(),
		     numBCDOF, topology::FieldBase::OTHER,
		     valueScale);
  } // if
  if (_dbRate) {
    const std::string& fieldLabel = 
      std::string("rate_") + std::string(fieldName);
    _parameters->add("rate", fieldLabel.c_str(),
		     numBCDOF, topology::FieldBase::OTHER,
		     rateScale);
    const std::string& timeLabel = 
      std::string("rate_time_") + std::string(fieldName);    
    _parameters->add("rate time", timeLabel.c_str(),
		     1, topology::FieldBase::SCALAR,
		     timeScale);
  } // if
  if (_dbChange) {
    const std::string& fieldLabel = 
      std::string("change_") + std::string(fieldName);
    _parameters->add("change", fieldLabel.c_str(),
		     numBCDOF, topology::FieldBase::OTHER,
		     valueScale);
    const std::string& timeLabel = 
      std::string("change_time_") + std::string(fieldName);
    _parameters->add("change time", timeLabel.c_str(),
		     1, topology::FieldBase::SCALAR,
		     timeScale);
  } // if
  _parameters->allocate(_points);
  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  
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
				 const double scale)
{ // _queryDB
  assert(name);
  assert(db);
  assert(_parameters);

  const topology::Mesh& mesh = _parameters->mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = _getNormalizer().lengthScale();

  double_array coordsVertex(spaceDim);
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  const int valueIndex = _parameters->sectionIndex(name);
  const int valueFiberDim = _parameters->sectionFiberDim(name);
  assert(valueIndex+valueFiberDim <= parametersFiberDim);
  double_array parametersVertex(parametersFiberDim);

  double_array valueVertex(querySize);

  const int numPoints = _points.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    // Get dimensionalized coordinates of vertex
    coordinates->restrictPoint(_points[iPoint], 
			       &coordsVertex[0], coordsVertex.size());
    _getNormalizer().dimensionalize(&coordsVertex[0], coordsVertex.size(),
				lengthScale);
    int err = db->query(&valueVertex[0], valueVertex.size(), 
			&coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Error querying for '" << name << "' at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") using spatial database " << db->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    _getNormalizer().nondimensionalize(&valueVertex[0], valueVertex.size(),
				   scale);

    // Update section
    assert(parametersFiberDim == 
	   parametersSection->getFiberDimension(_points[iPoint]));
    parametersSection->restrictPoint(_points[iPoint], &parametersVertex[0],
				     parametersVertex.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersVertex[valueIndex+i] = valueVertex[i];
    
    parametersSection->updatePoint(_points[iPoint], &parametersVertex[0]);
  } // for
} // _queryDB

// ----------------------------------------------------------------------
// Calculate temporal and spatial variation of value over the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValue(const double t)
{ // _calculateValue
  assert(_parameters);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const double timeScale = _getNormalizer().timeScale();

  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  double_array parametersVertex(parametersFiberDim);
  
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(numBCDOF == valueFiberDim);
  
  const int initialIndex = 
    (_dbInitial) ? _parameters->sectionIndex("initial") : -1;
  const int initialFiberDim = 
    (_dbInitial) ? _parameters->sectionFiberDim("initial") : 0;

  const int rateIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate") : -1;
  const int rateFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate") : 0;
  const int rateTimeIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate time") : -1;
  const int rateTimeFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate time") : 0;

  const int changeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change") : -1;
  const int changeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change") : 0;
  const int changeTimeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change time") : -1;
  const int changeTimeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change time") : 0;

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label
    
    assert(parametersFiberDim == parametersSection->getFiberDimension(p_bc));
    parametersSection->restrictPoint(p_bc, &parametersVertex[0],
				     parametersVertex.size());
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
      
      const double tRel = t - parametersVertex[rateTimeIndex];
      if (tRel > 0.0)  // rate of change integrated over time
	for (int iDim=0; iDim < numBCDOF; ++iDim)
	  parametersVertex[valueIndex+iDim] += 
	    parametersVertex[rateIndex+iDim] * tRel;
    } // if

    // Contribution from change of value
    if (_dbChange) {
      assert(changeIndex >= 0);
      assert(changeFiberDim == valueFiberDim);
      assert(changeTimeIndex >= 0);
      assert(changeTimeFiberDim == 1);

      const double tRel = t - parametersVertex[changeTimeIndex];
      if (tRel >= 0) { // change in value over time
	double scale = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = tRel;
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
	for (int iDim=0; iDim < numBCDOF; ++iDim)
	  parametersVertex[valueIndex+iDim] += 
	    parametersVertex[changeIndex+iDim] * scale;
      } // if
    } // if

    parametersSection->updatePoint(p_bc, &parametersVertex[0]);
  } // for
}  // _calculateValue

// ----------------------------------------------------------------------
// Calculate increment in temporal and spatial variation of value over
// the list of points.
void
pylith::bc::TimeDependentPoints::_calculateValueIncr(const double t0,
						     const double t1)
{ // _calculateValueIncr
  assert(_parameters);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  const double timeScale = _getNormalizer().timeScale();

  const ALE::Obj<RealUniformSection>& parametersSection = 
    _parameters->section();
  assert(!parametersSection.isNull());
  const int parametersFiberDim = _parameters->fiberDim();
  double_array parametersVertex(parametersFiberDim);
  
  const int valueIndex = _parameters->sectionIndex("value");
  const int valueFiberDim = _parameters->sectionFiberDim("value");
  assert(numBCDOF == valueFiberDim);
  
  const int rateIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate") : -1;
  const int rateFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate") : 0;
  const int rateTimeIndex = 
    (_dbRate) ? _parameters->sectionIndex("rate time") : -1;
  const int rateTimeFiberDim = 
    (_dbRate) ? _parameters->sectionFiberDim("rate time") : 0;

  const int changeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change") : -1;
  const int changeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change") : 0;
  const int changeTimeIndex = 
    (_dbChange) ? _parameters->sectionIndex("change time") : -1;
  const int changeTimeFiberDim = 
    (_dbChange) ? _parameters->sectionFiberDim("change time") : 0;

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const int p_bc = _points[iPoint]; // Get point label
    
    assert(parametersFiberDim == parametersSection->getFiberDimension(p_bc));
    parametersSection->restrictPoint(p_bc, &parametersVertex[0],
				     parametersVertex.size());
    for (int i=0; i < valueFiberDim; ++i)
      parametersVertex[valueIndex+i] = 0.0;

    // No contribution from initial value
    
    // Contribution from rate of change of value
    if (_dbRate) {
      assert(rateIndex >= 0);
      assert(rateFiberDim == valueFiberDim);
      assert(rateTimeIndex >= 0);
      assert(rateTimeFiberDim == 1);
      
      // Account for when rate dependence begins.
      const double tRate = parametersVertex[rateTimeIndex];
      double tIncr = 0.0;
      if (t0 > tRate) // rate dependence for t0 to t1
	tIncr = t1 - t0;
      else if (t1 > tRate) // rate dependence for tRef to t1
	tIncr = t1 - tRate;
      else
	tIncr = 0.0; // no rate dependence for t0 to t1

      if (tIncr > 0.0)  // rate of change integrated over time
	for (int iDim=0; iDim < numBCDOF; ++iDim)
	  parametersVertex[valueIndex+iDim] += 
	    parametersVertex[rateIndex+iDim] * tIncr;
    } // if
    
    // Contribution from change of value
    if (_dbChange) {
      assert(changeIndex >= 0);
      assert(changeFiberDim == valueFiberDim);
      assert(changeTimeIndex >= 0);
      assert(changeTimeFiberDim == 1);

      const double tChange = parametersVertex[changeTimeIndex];
      if (t0 >= tChange) { // increment is after change starts
	double scale0 = 1.0;
	double scale1 = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = t0 - tChange;
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
	for (int iDim=0; iDim < numBCDOF; ++iDim)
	  parametersVertex[valueIndex+iDim] += 
	    parametersVertex[changeIndex+iDim] * (scale1 - scale0);
      } else if (t1 >= tChange) { // increment spans when change starts
	double scale1 = 1.0;
	if (0 != _dbTimeHistory) {
	  double tDim = t1 - tChange;
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
	for (int iDim=0; iDim < numBCDOF; ++iDim)
	  parametersVertex[valueIndex+iDim] += 
	    parametersVertex[changeIndex+iDim] * scale1;
      } // if/else
    } // if

    parametersSection->updatePoint(p_bc, &parametersVertex[0]);
  } // for
}  // _calculateValueIncr


// End of file 
