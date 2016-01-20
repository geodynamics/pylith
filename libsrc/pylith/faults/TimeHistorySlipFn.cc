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

#include "TimeHistorySlipFn.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/faults/FaultCohesiveLagrange.hh" // USES isClampedVertex()

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

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
  PYLITH_METHOD_BEGIN;

  SlipTimeFn::deallocate();

  _dbAmplitude = 0; // :TODO: Use shared pointer
  _dbSlipTime = 0; // :TODO: Use shared pointer
  if (_dbTimeHistory)
    _dbTimeHistory->close();
  _dbTimeHistory = 0; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::TimeHistorySlipFn::initialize(const topology::Mesh& faultMesh,
					      const spatialdata::units::Nondimensional& normalizer,
					      const PylithScalar originTime)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_dbAmplitude);
  assert(_dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();

  // Get vertices in fault mesh
  PetscDM dmMesh = faultMesh.dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  delete _parameters; _parameters = new topology::Fields(faultMesh);assert(_parameters);

  _parameters->add("slip amplitude", "slip_amplitude");
  topology::Field& slipAmplitude = _parameters->get("slip amplitude");
  slipAmplitude.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slipAmplitude.allocate();
  slipAmplitude.scale(lengthScale);
  slipAmplitude.vectorFieldType(topology::FieldBase::VECTOR);
  topology::VecVisitorMesh slipAmplitudeVisitor(slipAmplitude);
  PetscScalar* slipAmplitudeArray = slipAmplitudeVisitor.localArray();

  _parameters->add("slip time", "slip_time");
  topology::Field& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipAmplitude, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

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
      std::ostringstream msg;
      msg << "Bad spatial dimension '" << spaceDim << " in TimeHistorySlipFn'." << std::endl;
      throw std::logic_error(msg.str());
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  // Get coordinates of vertices
  scalar_array vCoordsGlobal(spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar* coordsArray = coordsVisitor.localArray();

  PetscDM faultDMMesh = faultMesh.dmMesh();assert(faultDMMesh);
  PetscDMLabel clamped = NULL;
  PetscErrorCode err = DMGetLabel(faultDMMesh, "clamped", &clamped);PYLITH_CHECK_ERROR(err);

  _slipVertex.resize(spaceDim);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    if (FaultCohesiveLagrange::isClampedVertex(clamped, v)) {
      continue;
    } // if
      
    // Dimensionalize coordinates
    const PetscInt coff = coordsVisitor.sectionOffset(v);
    assert(spaceDim == coordsVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      vCoordsGlobal[d] = coordsArray[coff+d];
    } // for
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);

    // Slip amplitude
    const PetscInt saoff = slipAmplitudeVisitor.sectionOffset(v);
    assert(spaceDim == slipAmplitudeVisitor.sectionDof(v));        
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

    // Slip time
    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    assert(1 == slipTimeVisitor.sectionDof(v));
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

    for(PetscInt d = 0; d < spaceDim; ++d) {
      slipAmplitudeArray[saoff+d] = _slipVertex[d];
    } // for
    slipTimeArray[stoff] = _slipTimeVertex;
  } // for

  // Close databases.
  _dbAmplitude->close();
  _dbSlipTime->close();

  // Open time history database.
  _dbTimeHistory->open();
  _timeScale = timeScale;

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::TimeHistorySlipFn::slip(topology::Field* slip,
					const PylithScalar t)
{ // slip
  PYLITH_METHOD_BEGIN;

  assert(slip);
  assert(_parameters);
  assert(_dbTimeHistory);

  // Get vertices in fault mesh
  PetscDM dmMesh = _parameters->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Get sections
  const topology::Field& slipAmplitude = _parameters->get("slip amplitude");
  topology::VecVisitorMesh slipAmplitudeVisitor(slipAmplitude);
  const PetscScalar* slipAmplitudeArray = slipAmplitudeVisitor.localArray();

  const topology::Field& slipTime = _parameters->get("slip time");
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

  topology::VecVisitorMesh slipVisitor(*slip);
  PetscScalar* slipArray = slipVisitor.localArray();

  const int spaceDim = _slipVertex.size();
  PylithScalar amplitude = 0.0;
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt saoff = slipAmplitudeVisitor.sectionOffset(v);
    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    const PetscInt soff = slipVisitor.sectionOffset(v);

    assert(spaceDim == slipAmplitudeVisitor.sectionDof(v));
    assert(1 == slipTimeVisitor.sectionDof(v));
    assert(spaceDim == slipVisitor.sectionDof(v));

    PylithScalar relTime = t - slipTimeArray[stoff];
    if (relTime >= 0.0) {
      relTime *= _timeScale;
      const int err = _dbTimeHistory->query(&amplitude, relTime);
      if (err) {
        std::ostringstream msg;
        msg << "Error querying for time '" << relTime
            << "' in time history database '"
            << _dbTimeHistory->label() << "'.";
        throw std::runtime_error(msg.str());
      } // if

      for(PetscInt d = 0; d < spaceDim; ++d) {
        slipArray[soff+d] += slipAmplitudeArray[saoff+d] * amplitude;
      } // for
    } // else
  } // for

  PetscLogFlops((vEnd-vStart) * 3);

  PYLITH_METHOD_END;
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field&
pylith::faults::TimeHistorySlipFn::finalSlip(void)
{ // finalSlip
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_parameters->get("slip amplitude"));
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field&
pylith::faults::TimeHistorySlipFn::slipTime(void)
{ // slipTime
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_parameters->get("slip time"));
} // slipTime


// End of file 
