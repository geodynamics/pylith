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

#include "BruneSlipFn.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/faults/FaultCohesiveLagrange.hh" // USES isClampedVertex()

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::BruneSlipFn::BruneSlipFn(void) :
  _slipTimeVertex(0),
  _riseTimeVertex(0),
  _dbFinalSlip(0),
  _dbSlipTime(0),
  _dbRiseTime(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::BruneSlipFn::~BruneSlipFn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::BruneSlipFn::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  SlipTimeFn::deallocate();

  _dbFinalSlip = 0; // :TODO: Use shared pointer.
  _dbSlipTime = 0; // :TODO: Use shared pointer.
  _dbRiseTime = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::BruneSlipFn::initialize(const topology::Mesh& faultMesh,
					const spatialdata::units::Nondimensional& normalizer,
					const PylithScalar originTime)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_dbFinalSlip);
  assert(_dbSlipTime);
  assert(_dbRiseTime);

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
  _parameters->add("final slip", "final_slip");
  topology::Field& finalSlip = _parameters->get("final slip");
  finalSlip.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  finalSlip.allocate();
  finalSlip.scale(lengthScale);
  finalSlip.vectorFieldType(topology::FieldBase::VECTOR);
  topology::VecVisitorMesh finalSlipVisitor(finalSlip);
  PetscScalar* finalSlipArray = finalSlipVisitor.localArray();

  _parameters->add("slip time", "slip_time");
  topology::Field& slipTime = _parameters->get("slip time");
  slipTime.newSection(finalSlip, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

  _parameters->add("rise time", "rise_time");
  topology::Field& riseTime = _parameters->get("rise time");
  riseTime.newSection(finalSlip, 1);
  riseTime.allocate();
  riseTime.scale(timeScale);
  riseTime.vectorFieldType(topology::FieldBase::SCALAR);
  topology::VecVisitorMesh riseTimeVisitor(riseTime);
  PetscScalar* riseTimeArray = riseTimeVisitor.localArray();

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
      std::ostringstream msg;
      msg << "Bad spatial dimension '" << spaceDim << "' in BruneSlipFn." << std::endl;
      throw std::logic_error(msg.str());
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  _dbRiseTime->open();
  const char* riseTimeValues[] = {"rise-time"};
  _dbRiseTime->queryVals(riseTimeValues, 1);

  // Get coordinates of vertices
  scalar_array vCoordsGlobal(spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar* coordsArray = coordsVisitor.localArray();

  PetscDM faultDMMesh = faultMesh.dmMesh();assert(faultDMMesh);
  PetscDMLabel clamped = NULL;
  PetscErrorCode err = DMGetLabel(faultDMMesh, "clamped", &clamped);PYLITH_CHECK_ERROR(err);

  _slipVertex.resize(spaceDim);
  for (PetscInt v = vStart; v < vEnd; ++v) {
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
    
    // Final slip
    const PetscInt fsoff = finalSlipVisitor.sectionOffset(v);
    assert(spaceDim == finalSlipVisitor.sectionDof(v));
    err = _dbFinalSlip->query(&_slipVertex[0], _slipVertex.size(), &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database '" << _dbFinalSlip->label() << "'.";
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
      msg << ") using spatial database '" << _dbSlipTime->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipTimeVertex, 1, timeScale);
    // add origin time to rupture time
    _slipTimeVertex += originTime;

    // Rise time
    const PetscInt rtoff = riseTimeVisitor.sectionOffset(v);
    assert(1 == riseTimeVisitor.sectionDof(v));
    err = _dbRiseTime->query(&_riseTimeVertex, 1, &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find rise time at (";
      for (int i=0; i < spaceDim; ++i)
        msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database '" << _dbRiseTime->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_riseTimeVertex, 1, timeScale);

    for(PetscInt d = 0; d < spaceDim; ++d) {
      finalSlipArray[fsoff+d] = _slipVertex[d];
    }
    slipTimeArray[stoff] = _slipTimeVertex;
    riseTimeArray[rtoff] = _riseTimeVertex;
  } // for

  // Close databases
  _dbFinalSlip->close();
  _dbSlipTime->close();
  _dbRiseTime->close();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::BruneSlipFn::slip(topology::Field* slip,
				  const PylithScalar t)
{ // slip
  PYLITH_METHOD_BEGIN;

  assert(slip);
  assert(_parameters);

  // Get vertices in fault mesh
  PetscDM dmMesh = _parameters->mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  // Get sections
  const topology::Field& finalSlip = _parameters->get("final slip");
  topology::VecVisitorMesh finalSlipVisitor(finalSlip);
  const PetscScalar* finalSlipArray = finalSlipVisitor.localArray();

  const topology::Field& slipTime = _parameters->get("slip time");
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

  const topology::Field& riseTime = _parameters->get("rise time");
  topology::VecVisitorMesh riseTimeVisitor(riseTime);
  const PetscScalar* riseTimeArray = riseTimeVisitor.localArray();

  topology::VecVisitorMesh slipVisitor(*slip);
  PetscScalar* slipArray = slipVisitor.localArray();
  DMLabel clamped = NULL;
  PetscErrorCode err;

  err = DMGetLabel(dmMesh, "clamped", &clamped);PYLITH_CHECK_ERROR(err);
  const int spaceDim = _slipVertex.size();
  for (PetscInt v = vStart; v < vEnd; ++v) {
    if (FaultCohesiveLagrange::isClampedVertex(clamped, v)) {
      continue;
    } // if

    const PetscInt fsoff = finalSlipVisitor.sectionOffset(v);
    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    const PetscInt rtoff = riseTimeVisitor.sectionOffset(v);
    const PetscInt soff = slipVisitor.sectionOffset(v);

    assert(spaceDim == finalSlipVisitor.sectionDof(v));
    assert(1 == slipTimeVisitor.sectionDof(v));
    assert(1 == riseTimeVisitor.sectionDof(v));
    assert(spaceDim == slipVisitor.sectionDof(v));

    PylithScalar finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += finalSlipArray[fsoff+i]*finalSlipArray[fsoff+i];
    finalSlipMag = sqrt(finalSlipMag);

    const PylithScalar slip = _slipFn(t-slipTimeArray[stoff], finalSlipMag, riseTimeArray[rtoff]);
    const PylithScalar scale = finalSlipMag > 0.0 ? slip / finalSlipMag : 0.0;
    
    // Update field
    for(PetscInt d = 0; d < spaceDim; ++d) {
      slipArray[soff+d] += finalSlipArray[fsoff+d] * scale;
    } // for
  } // for

  PetscLogFlops((vEnd-vStart) * (2+8 + 3*spaceDim));

  PYLITH_METHOD_END;
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field&
pylith::faults::BruneSlipFn::finalSlip(void)
{ // finalSlip
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_parameters->get("final slip"));
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field&
pylith::faults::BruneSlipFn::slipTime(void)
{ // slipTime
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_parameters->get("slip time"));
} // slipTime


// End of file 
