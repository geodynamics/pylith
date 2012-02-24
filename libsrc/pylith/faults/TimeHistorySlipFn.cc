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

#include "TimeHistorySlipFn.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh::label_sequence label_sequence;
typedef pylith::topology::SubMesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::TimeHistorySlipFn::TimeHistorySlipFn(void) :
  _slipTimeVertex(0),
  _timeScale(1.0),
  _dbAmplitude(0),
  _dbSlipTime(0),
  _dbTimeHistory(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::TimeHistorySlipFn::~TimeHistorySlipFn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::TimeHistorySlipFn::deallocate(void)
{ // deallocate
  SlipTimeFn::deallocate();

  _dbAmplitude = 0; // :TODO: Use shared pointer
  _dbSlipTime = 0; // :TODO: Use shared pointer
  if (0 != _dbTimeHistory)
    _dbTimeHistory->close();
  _dbTimeHistory = 0; // :TODO: Use shared pointer
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::TimeHistorySlipFn::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer,
			    const double originTime)
{ // initialize
  assert(0 != _dbAmplitude);
  assert(0 != _dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = normalizer.lengthScale();
  const double timeScale = normalizer.timeScale();

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = faultMesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(0 != _parameters);
  _parameters->add("slip amplitude", "slip_amplitude");

  topology::Field<topology::SubMesh>& slipAmplitude = 
    _parameters->get("slip amplitude");
  slipAmplitude.newSection(vertices, spaceDim);
  slipAmplitude.allocate();
  slipAmplitude.scale(lengthScale);
  slipAmplitude.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& slipAmplitudeSection = slipAmplitude.section();
  assert(!slipAmplitudeSection.isNull());  

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipAmplitude, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());

  logger.stagePop();

  // Open databases and set query values
  _dbAmplitude->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"fault-opening"};
      _dbAmplitude->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"left-lateral-slip", "fault-opening"};
      _dbAmplitude->queryVals(slipValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* slipValues[] = {"left-lateral-slip", "reverse-slip", 
				  "fault-opening"};
      _dbAmplitude->queryVals(slipValues, 3);
      break;
    } // case 3
    default :
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in TimeHistorySlipFn.");
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  // Get coordinates of vertices
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  _slipVertex.resize(spaceDim);
  double_array vCoordsGlobal(spaceDim);  
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
			      lengthScale);
        
    int err = _dbAmplitude->query(&_slipVertex[0], _slipVertex.size(), 
				  &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip amplitude at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbAmplitude->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipVertex[0], _slipVertex.size(),
				 lengthScale);

    err = _dbSlipTime->query(&_slipTimeVertex, 1, 
			     &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip initiation time at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbSlipTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipTimeVertex, 1, timeScale);
    // add origin time to rupture time
    _slipTimeVertex += originTime;

    slipAmplitudeSection->updatePoint(*v_iter, &_slipVertex[0]);
    slipTimeSection->updatePoint(*v_iter, &_slipTimeVertex);
  } // for

  // Close databases.
  _dbAmplitude->close();
  _dbSlipTime->close();

  // Open time history database.
  _dbTimeHistory->open();
  _timeScale = timeScale;
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::TimeHistorySlipFn::slip(topology::Field<topology::SubMesh>* slip,
				 const double t)
{ // slip
  assert(0 != slip);
  assert(0 != _parameters);
  assert(0 != _dbTimeHistory);

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = slip->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const topology::Field<topology::SubMesh>& slipAmplitude = 
    _parameters->get("slip amplitude");
  const ALE::Obj<RealSection>& slipAmplitudeSection = slipAmplitude.section();
  assert(!slipAmplitudeSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  double amplitude = 0.0;
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    slipAmplitudeSection->restrictPoint(*v_iter, 
					&_slipVertex[0], _slipVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);

    double relTime = t - _slipTimeVertex;
    if (relTime < 0.0)
      _slipVertex = 0.0;
    else {
      relTime *= _timeScale;
      const int err = _dbTimeHistory->query(&amplitude, relTime);
      if (0 != err) {
	std::ostringstream msg;
	msg << "Error querying for time '" << relTime
	    << "' in time history database "
	    << _dbTimeHistory->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      _slipVertex *= amplitude;
    } // else
    
    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * 3);
} // slip

// ----------------------------------------------------------------------
// Get increment of slip on fault surface between time t0 and t1.
void
pylith::faults::TimeHistorySlipFn::slipIncr(topology::Field<topology::SubMesh>* slip,
				     const double t0,
				     const double t1)
{ // slipIncr
  assert(0 != slip);
  assert(0 != _parameters);
  assert(0 != _dbTimeHistory);

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = slip->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const topology::Field<topology::SubMesh>& slipAmplitude = 
    _parameters->get("slip amplitude");
  const ALE::Obj<RealSection>& slipAmplitudeSection = slipAmplitude.section();
  assert(!slipAmplitudeSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  double amplitude0 = 0.0;
  double amplitude1 = 0.0;
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    slipAmplitudeSection->restrictPoint(*v_iter, &_slipVertex[0], _slipVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);

    double relTime0 = t0 - _slipTimeVertex;
    double relTime1 = t1 - _slipTimeVertex;
    if (relTime1 < 0.0)
      _slipVertex = 0.0;
    else {
      relTime0 *= _timeScale;
      relTime1 *= _timeScale;
      int err = _dbTimeHistory->query(&amplitude0, relTime0);
      if (0 != err) {
	std::ostringstream msg;
	msg << "Error querying for time '" << relTime0
	    << "' in time history database "
	    << _dbTimeHistory->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      err = _dbTimeHistory->query(&amplitude1, relTime1);
      if (0 != err) {
	std::ostringstream msg;
	msg << "Error querying for time '" << relTime1
	    << "' in time history database "
	    << _dbTimeHistory->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      _slipVertex *= amplitude1 - amplitude0;
    } // else

    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * 6);
} // slipIncr

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::TimeHistorySlipFn::finalSlip(void)
{ // finalSlip
  return _parameters->get("slip amplitude");
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::TimeHistorySlipFn::slipTime(void)
{ // slipTime
  return _parameters->get("slip time");
} // slipTime


// End of file 
