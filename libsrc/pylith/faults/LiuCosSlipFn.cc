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

#include "LiuCosSlipFn.hh" // implementation of object methods

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
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh::label_sequence label_sequence;
typedef pylith::topology::SubMesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::LiuCosSlipFn::LiuCosSlipFn(void) :
  _slipTimeVertex(0),
  _riseTimeVertex(0),
  _dbFinalSlip(0),
  _dbSlipTime(0),
  _dbRiseTime(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::LiuCosSlipFn::~LiuCosSlipFn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::LiuCosSlipFn::deallocate(void)
{ // deallocate
  SlipTimeFn::deallocate();

  _dbFinalSlip = 0; // :TODO: Use shared pointer.
  _dbSlipTime = 0; // :TODO: Use shared pointer.
  _dbRiseTime = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::LiuCosSlipFn::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer,
			    const PylithScalar originTime)
{ // initialize
  assert(0 != _dbFinalSlip);
  assert(0 != _dbSlipTime);
  assert(0 != _dbRiseTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();

  // Get vertices in fault mesh
  DM             dmMesh = faultMesh.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(0 != _parameters);
  _parameters->add("final slip", "final_slip");
  topology::Field<topology::SubMesh>& finalSlip = _parameters->get("final slip");
  finalSlip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  finalSlip.allocate();
  finalSlip.scale(lengthScale);
  finalSlip.vectorFieldType(topology::FieldBase::VECTOR);
  PetscSection finalSlipSection  = finalSlip.petscSection();
  Vec          finalSlipVec      = finalSlip.localVector();
  PetscScalar *finalSlipArray;
  assert(finalSlipSection);assert(finalSlipVec);

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(finalSlip, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  PetscSection slipTimeSection  = slipTime.petscSection();
  Vec          slipTimeVec      = slipTime.localVector();
  PetscScalar *slipTimeArray;
  assert(slipTimeSection);assert(slipTimeVec);

  _parameters->add("rise time", "rise_time");
  topology::Field<topology::SubMesh>& riseTime = _parameters->get("rise time");
  riseTime.newSection(finalSlip, 1);
  riseTime.allocate();
  riseTime.scale(timeScale);
  riseTime.vectorFieldType(topology::FieldBase::SCALAR);
  PetscSection riseTimeSection  = riseTime.petscSection();
  Vec          riseTimeVec      = riseTime.localVector();
  PetscScalar *riseTimeArray;
  assert(riseTimeSection);assert(riseTimeVec);

  logger.stagePop();

  // Open databases and set query values
  _dbFinalSlip->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* slipValues[] = {"fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 1);
      break;
    } // case 1
    case 2 : {
      const char* slipValues[] = {"left-lateral-slip", "fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 2);
      break;
    } // case 2
    case 3 : {
      const char* slipValues[] = {"left-lateral-slip", "reverse-slip", 
				  "fault-opening"};
      _dbFinalSlip->queryVals(slipValues, 3);
      break;
    } // case 3
    default :
      std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad spatial dimension in LiuCosSlipFn.");
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  _dbRiseTime->open();
  const char* riseTimeValues[] = {"rise-time"};
  _dbRiseTime->queryVals(riseTimeValues, 1);

  // Get coordinates of vertices
  scalar_array vCoordsGlobal(spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coordArray;
  err = DMComplexGetCoordinateSection(dmMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(dmMesh, &coordVec);CHECK_PETSC_ERROR(err);
  assert(coordSection);assert(coordVec);

  _slipVertex.resize(spaceDim);
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec,  &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(riseTimeVec,  &riseTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt fsdof, fsoff, stdof, stoff, rtdof, rtoff, cdof, coff;

    err = PetscSectionGetDof(finalSlipSection, v, &fsdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &fsoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(riseTimeSection, v, &rtdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(riseTimeSection, v, &rtoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(coordSection, v, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v, &coff);CHECK_PETSC_ERROR(err);
    assert(cdof == spaceDim);

    for(PetscInt d = 0; d < cdof; ++d) {
      vCoordsGlobal[d] = coordArray[coff+d];
    }
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);
    
    int err = _dbFinalSlip->query(&_slipVertex[0], _slipVertex.size(), &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbFinalSlip->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipVertex[0], _slipVertex.size(), lengthScale);

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

    err = _dbRiseTime->query(&_riseTimeVertex, 1, &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find rise time at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbRiseTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_riseTimeVertex, 1, timeScale);

    for(PetscInt d = 0; d < fsdof; ++d) {
      finalSlipArray[fsoff+d] = _slipVertex[d];
    }
    slipTimeArray[stoff] = _slipTimeVertex;
    riseTimeArray[stoff] = _riseTimeVertex;
  } // for
  err = VecRestoreArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(slipTimeVec,  &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(riseTimeVec,  &riseTimeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(coordVec, &coordArray);CHECK_PETSC_ERROR(err);

  // Close databases
  _dbFinalSlip->close();
  _dbSlipTime->close();
  _dbRiseTime->close();
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::LiuCosSlipFn::slip(topology::Field<topology::SubMesh>* slip,
				  const PylithScalar t)
{ // slip
  assert(0 != slip);
  assert(0 != _parameters);

  // Get vertices in fault mesh
  DM             dmMesh = slip->mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  // Get sections
  const topology::Field<topology::SubMesh>& finalSlip = _parameters->get("final slip");
  PetscSection finalSlipSection  = finalSlip.petscSection();
  Vec          finalSlipVec      = finalSlip.localVector();
  PetscScalar *finalSlipArray;
  assert(finalSlipSection);assert(finalSlipVec);
  const topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  PetscSection slipTimeSection  = slipTime.petscSection();
  Vec          slipTimeVec      = slipTime.localVector();
  PetscScalar *slipTimeArray;
  assert(slipTimeSection);assert(slipTimeVec);
  const topology::Field<topology::SubMesh>& riseTime = _parameters->get("rise time");
  PetscSection riseTimeSection  = riseTime.petscSection();
  Vec          riseTimeVec      = riseTime.localVector();
  PetscScalar *riseTimeArray;
  assert(riseTimeSection);assert(riseTimeVec);
  PetscSection slipSection  = slip->petscSection();
  Vec          slipVec      = slip->localVector();
  PetscScalar *slipArray;
  assert(slipSection);assert(slipVec);

  const int spaceDim = _slipVertex.size();
  err = VecGetArray(finalSlipVec, &finalSlipArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipTimeVec, &slipTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(riseTimeVec, &riseTimeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(slipVec, &slipArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt fsdof, fsoff, stdof, stoff, rtdof, rtoff, sdof, soff;

    err = PetscSectionGetDof(finalSlipSection, v, &fsdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(finalSlipSection, v, &fsoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipTimeSection, v, &stdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipTimeSection, v, &stoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(riseTimeSection, v, &rtdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(riseTimeSection, v, &rtoff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(slipSection, v, &sdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(slipSection, v, &soff);CHECK_PETSC_ERROR(err);
    assert(fsdof == sdof);
    assert(stdof == 1);
    assert(rtdof == 1);

    PylithScalar finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += finalSlipArray[fsoff+i]*finalSlipArray[fsoff+i];
    finalSlipMag = sqrt(finalSlipMag);

    const PylithScalar slip = _slipFn(t-slipTimeArray[stoff], finalSlipMag, riseTimeArray[rtoff]);
    const PylithScalar scale = finalSlipMag > 0.0 ? slip / finalSlipMag : 0.0;
    
    // Update field
    for(PetscInt d = 0; d < sdof; ++d) {
      slipArray[soff+d] += finalSlipArray[fsoff+d] * scale;
    }
  } // for

  PetscLogFlops((vEnd-vStart) * (2+28 + 3*_slipVertex.size()));
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::LiuCosSlipFn::finalSlip(void)
{ // finalSlip
  return _parameters->get("final slip");
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::LiuCosSlipFn::slipTime(void)
{ // slipTime
  return _parameters->get("slip time");
} // slipTime


// End of file 
