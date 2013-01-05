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
			    const PylithScalar originTime)
{ // initialize
  assert(0 != _dbAmplitude);
  assert(0 != _dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();

  // Memory logging
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  // Get vertices in fault mesh
  DM             dmMesh = faultMesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(0 != _parameters);
  _parameters->add("slip amplitude", "slip_amplitude");

  topology::Field<topology::SubMesh>& slipAmplitude = _parameters->get("slip amplitude");
  slipAmplitude.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slipAmplitude.allocate();
  slipAmplitude.scale(lengthScale);
  slipAmplitude.vectorFieldType(topology::FieldBase::VECTOR);
  PetscSection finalSlipSection  = slipAmplitude.petscSection();
  Vec          finalSlipVec      = slipAmplitude.localVector();
  PetscScalar *finalSlipArray;
  assert(finalSlipSection);assert(finalSlipVec);

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipAmplitude, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  PetscSection slipTimeSection  = slipTime.petscSection();
  Vec          slipTimeVec      = slipTime.localVector();
  PetscScalar *slipTimeArray;
  assert(slipTimeSection);assert(slipTimeVec);

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
  scalar_array vCoordsGlobal(spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coordArray;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  _slipVertex.resize(spaceDim);
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec,  &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt fsdof, fsoff, stdof, stoff, cdof, coff;

    err = PetscSectionGetDof(finalSlipSection, v, &fsdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &fsoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
    assert(cdof == spaceDim);

    for(PetscInt d = 0; d < cdof; ++d) {
      vCoordsGlobal[d] = coordArray[coff+d];
    }
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);
        
    int err = _dbAmplitude->query(&_slipVertex[0], _slipVertex.size(), &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip amplitude at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbAmplitude->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipVertex[0], _slipVertex.size(), lengthScale);

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

    for(PetscInt d = 0; d < fsdof; ++d) {
      finalSlipArray[fsoff+d] = _slipVertex[d];
    }
    slipTimeArray[stoff] = _slipTimeVertex;
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
				 const PylithScalar t)
{ // slip
  assert(0 != slip);
  assert(0 != _parameters);
  assert(0 != _dbTimeHistory);

  // Get vertices in fault mesh
  DM             dmMesh = slip->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Get sections
  const topology::Field<topology::SubMesh>& finalSlip = _parameters->get("slip amplitude");
  PetscSection finalSlipSection  = finalSlip.petscSection();
  Vec          finalSlipVec      = finalSlip.localVector();
  PetscScalar *finalSlipArray;
  assert(finalSlipSection);assert(finalSlipVec);
  const topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  PetscSection slipTimeSection  = slipTime.petscSection();
  Vec          slipTimeVec      = slipTime.localVector();
  PetscScalar *slipTimeArray;
  assert(slipTimeSection);assert(slipTimeVec);
  PetscSection slipSection  = slip->petscSection();
  Vec          slipVec      = slip->localVector();
  PetscScalar *slipArray;
  assert(slipSection);assert(slipVec);

  PylithScalar amplitude = 0.0;
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt fsdof, fsoff, stdof, stoff, sdof, soff;

    err = PetscSectionGetDof(finalSlipSection, v, &fsdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &fsoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipSection, v, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &soff);CHECK_PETSC_ERROR(err);
    assert(fsdof == sdof);
    assert(stdof == 1);

    PylithScalar relTime = t - slipTimeArray[stoff];
    if (relTime >= 0.0) {
      relTime *= _timeScale;
      const int err = _dbTimeHistory->query(&amplitude, relTime);
      if (0 != err) {
        std::ostringstream msg;
        msg << "Error querying for time '" << relTime
            << "' in time history database "
            << _dbTimeHistory->label() << ".";
        throw std::runtime_error(msg.str());
      } // if

      for(PetscInt d = 0; d < sdof; ++d) {
        slipArray[soff+d] += finalSlipArray[fsoff+d] * amplitude;
      }
    } // else
  } // for

  PetscLogFlops((vEnd-vStart) * 3);
} // slip

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
