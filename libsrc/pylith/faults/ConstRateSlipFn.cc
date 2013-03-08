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

#include "ConstRateSlipFn.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::ConstRateSlipFn::ConstRateSlipFn(void) :
  _slipTimeVertex(0),
  _dbSlipRate(0),
  _dbSlipTime(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::ConstRateSlipFn::~ConstRateSlipFn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::ConstRateSlipFn::deallocate(void)
{ // deallocate
  SlipTimeFn::deallocate();

  _dbSlipRate = 0; // :TODO: Use shared pointer.
  _dbSlipTime = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::ConstRateSlipFn::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer,
			    const PylithScalar originTime)
{ // initialize
  assert(_dbSlipRate);
  assert(_dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();
  const PylithScalar velocityScale =
    normalizer.lengthScale() / normalizer.timeScale();

  // Get vertices in fault mesh
  PetscDM dmMesh = faultMesh.dmMesh();assert(dmMesh);
  PetscInt vStart, vEnd;
  PetscErrorCode err;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(_parameters);
  _parameters->add("slip rate", "slip_rate");
  topology::Field<topology::SubMesh>& slipRate = _parameters->get("slip rate");
  slipRate.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slipRate.allocate();
  slipRate.scale(velocityScale);
  slipRate.vectorFieldType(topology::FieldBase::VECTOR);
  PetscSection slipRateSection = slipRate.petscSection();assert(slipRateSection);
  PetscVec slipRateVec = slipRate.localVector();assert(slipRateVec);
  PetscScalar *slipRateArray = NULL;

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipRate, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  PetscSection slipTimeSection = slipTime.petscSection();assert(slipTimeSection);
  PetscVec slipTimeVec = slipTime.localVector();assert(slipTimeVec);
  PetscScalar *slipTimeArray = NULL;

  logger.stagePop();

  // Open databases and set query values
  _dbSlipRate->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"fault-opening"};
      _dbSlipRate->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"left-lateral-slip", "fault-opening"};
      _dbSlipRate->queryVals(slipValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* slipValues[] = {"left-lateral-slip", "reverse-slip", 
				  "fault-opening"};
      _dbSlipRate->queryVals(slipValues, 3);
      break;
    } // case 3
    default :
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in ConstRateSlipFn.");
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  // Get coordinates of vertices
  scalar_array vCoordsGlobal(spaceDim);
  PetscSection coordSection = NULL;
  PetscVec coordVec = NULL;
  PetscScalar *coordArray = NULL;
  err = DMPlexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);assert(coordSection);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);assert(coordVec);

  _slipRateVertex.resize(spaceDim);
  err = VecGetArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt srdof, sroff, stdof, stoff, cdof, coff;
    err = PetscSectionGetDof(slipRateSection, v, &srdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipRateSection, v, &sroff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
    assert(cdof == spaceDim);

    for(PetscInt d = 0; d < cdof; ++d) {
      vCoordsGlobal[d] = coordArray[coff+d];
    } // for
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);
    
    int err = _dbSlipRate->query(&_slipRateVertex[0], _slipRateVertex.size(), &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbSlipRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipRateVertex[0], _slipRateVertex.size(), velocityScale);

    err = _dbSlipTime->query(&_slipTimeVertex, 1, &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
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

    for(PetscInt d = 0; d < srdof; ++d) {
      slipRateArray[sroff+d] = _slipRateVertex[d];
    } // for
    slipTimeArray[stoff] = _slipTimeVertex;
  } // for
  err = VecRestoreArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);

  // Close databases
  _dbSlipRate->close();
  _dbSlipTime->close();
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::ConstRateSlipFn::slip(topology::Field<topology::SubMesh>* slip,
				      const PylithScalar t)
{ // slip
  assert(slip);
  assert(_parameters);

  // Get vertices in fault mesh
  PetscDM dmMesh = slip->mesh().dmMesh();assert(dmMesh);
  PetscInt vStart, vEnd;
  PetscErrorCode err;
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Get sections
  const topology::Field<topology::SubMesh>& slipRate = _parameters->get("slip rate");
  PetscSection slipRateSection = slipRate.petscSection();assert(slipRateSection);
  PetscVec slipRateVec = slipRate.localVector();assert(slipRateVec);
  PetscScalar *slipRateArray = NULL;
  err = VecGetArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);

  const topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  PetscSection slipTimeSection = slipTime.petscSection();assert(slipTimeSection);
  PetscVec slipTimeVec = slipTime.localVector();assert(slipTimeVec);
  PetscScalar *slipTimeArray = NULL;
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);

  PetscSection slipSection = slip->petscSection();assert(slipSection);
  PetscVec slipVec = slip->localVector();assert(slipVec);
  PetscScalar *slipArray = NULL;
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt srdof, sroff, stdof, stoff, sdof, soff;
    err = PetscSectionGetDof(slipRateSection, v, &srdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipRateSection, v, &sroff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipSection, v, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &soff);CHECK_PETSC_ERROR(err);
    assert(srdof == sdof);
    assert(stdof == 1);

    const PylithScalar relTime = t - slipTimeArray[stoff];
    if (relTime > 0.0) {
      for(PetscInt d = 0; d < sdof; ++d) {
        slipArray[soff+d] += slipRateArray[sroff+d] * relTime; // Convert slip rate to slip
      } // for
    } // if
  } // for
  err = VecRestoreArray(slipRateVec, &slipRateArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);

  PetscLogFlops((vEnd-vStart) * (1 + _slipRateVertex.size()));
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::ConstRateSlipFn::finalSlip(void)
{ // finalSlip
  // Slip rate is parameter instead of final slip.
  return _parameters->get("slip rate");
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::ConstRateSlipFn::slipTime(void)
{ // slipTime
  return _parameters->get("slip time");
} // slipTime


// End of file 
