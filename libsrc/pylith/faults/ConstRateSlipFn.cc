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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "ConstRateSlipFn.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum

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
  PYLITH_METHOD_BEGIN;

  SlipTimeFn::deallocate();

  _dbSlipRate = 0; // :TODO: Use shared pointer.
  _dbSlipTime = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::ConstRateSlipFn::initialize(const topology::Mesh& faultMesh,
					    const spatialdata::units::Nondimensional& normalizer,
					    const PylithScalar originTime)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(_dbSlipRate);
  assert(_dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  const PylithScalar lengthScale = normalizer.lengthScale();
  const PylithScalar timeScale = normalizer.timeScale();
  const PylithScalar velocityScale = normalizer.lengthScale() / normalizer.timeScale();

  // Get vertices in fault mesh
  PetscDM dmMesh = faultMesh.dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  delete _parameters; _parameters = new topology::Fields(faultMesh);assert(_parameters);

  _parameters->add("slip rate", "slip_rate");
  topology::Field& slipRate = _parameters->get("slip rate");
  slipRate.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  slipRate.allocate();
  slipRate.scale(velocityScale);
  slipRate.vectorFieldType(topology::FieldBase::VECTOR);
  topology::VecVisitorMesh slipRateVisitor(slipRate);
  PetscScalar* slipRateArray = slipRateVisitor.localArray();

  _parameters->add("slip time", "slip_time");
  topology::Field& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipRate, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

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
      std::ostringstream msg;
      msg << "Bad spatial dimension '" << spaceDim << " in ConstRateSlipFn'." << std::endl;
      throw std::logic_error(msg.str());
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  // Get coordinates of vertices
  scalar_array vCoordsGlobal(spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);
  PetscScalar* coordsArray = coordsVisitor.localArray();

  _slipRateVertex.resize(spaceDim);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    // Dimensionalize coordinates
    const PetscInt coff = coordsVisitor.sectionOffset(v);
    assert(spaceDim == coordsVisitor.sectionDof(v));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      vCoordsGlobal[d] = coordsArray[coff+d];
    } // for
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(), lengthScale);
    
    // Slip Rate
    const PetscInt sroff = slipRateVisitor.sectionOffset(v);
    assert(spaceDim == slipRateVisitor.sectionDof(v));
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
      slipRateArray[sroff+d] = _slipRateVertex[d];
    } // for
    slipTimeArray[stoff] = _slipTimeVertex;
  } // for

  // Close databases
  _dbSlipRate->close();
  _dbSlipTime->close();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::ConstRateSlipFn::slip(topology::Field* slip,
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
  const topology::Field& slipRate = _parameters->get("slip rate");
  topology::VecVisitorMesh slipRateVisitor(slipRate);
  const PetscScalar* slipRateArray = slipRateVisitor.localArray();

  const topology::Field& slipTime = _parameters->get("slip time");
  topology::VecVisitorMesh slipTimeVisitor(slipTime);
  const PetscScalar* slipTimeArray = slipTimeVisitor.localArray();

  topology::VecVisitorMesh slipVisitor(*slip);
  PetscScalar* slipArray = slipVisitor.localArray();

  const int spaceDim = _slipRateVertex.size();
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt sroff = slipRateVisitor.sectionOffset(v);
    const PetscInt stoff = slipTimeVisitor.sectionOffset(v);
    const PetscInt soff = slipVisitor.sectionOffset(v);

    assert(spaceDim == slipRateVisitor.sectionDof(v));
    assert(1 == slipTimeVisitor.sectionDof(v));
    assert(spaceDim == slipVisitor.sectionDof(v));

    const PylithScalar relTime = t - slipTimeArray[stoff];
    if (relTime > 0.0) {
      for(PetscInt d = 0; d < spaceDim; ++d) {
        slipArray[soff+d] += slipRateArray[sroff+d] * relTime; // Convert slip rate to slip
      } // for
    } // if
  } // for

  PetscLogFlops((vEnd-vStart) * (1 + _slipRateVertex.size()));

  PYLITH_METHOD_END;
} // slip

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field&
pylith::faults::ConstRateSlipFn::finalSlip(void)
{ // finalSlip
  PYLITH_METHOD_BEGIN;

  // Slip rate is parameter instead of final slip.
  PYLITH_METHOD_RETURN(_parameters->get("slip rate"));
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field&
pylith::faults::ConstRateSlipFn::slipTime(void)
{ // slipTime
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_RETURN(_parameters->get("slip time"));
} // slipTime


// End of file 
