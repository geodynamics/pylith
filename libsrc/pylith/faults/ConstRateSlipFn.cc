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
// Copyright (c) 2010-2011 University of California, Davis
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
typedef pylith::topology::SubMesh::SieveMesh SieveMesh;
typedef pylith::topology::SubMesh::SieveMesh::label_sequence label_sequence;
typedef pylith::topology::SubMesh::RealSection RealSection;

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
			    const double originTime)
{ // initialize
  assert(0 != _dbSlipRate);
  assert(0 != _dbSlipTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = normalizer.lengthScale();
  const double timeScale = normalizer.timeScale();
  const double velocityScale =
    normalizer.lengthScale() / normalizer.timeScale();

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = faultMesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Fault");

  delete _parameters; _parameters = new topology::Fields<topology::Field<topology::SubMesh> >(faultMesh);
  assert(0 != _parameters);
  _parameters->add("slip rate", "slip_rate");
  topology::Field<topology::SubMesh>& slipRate = _parameters->get("slip rate");
  slipRate.newSection(vertices, spaceDim);
  slipRate.allocate();
  slipRate.scale(velocityScale);
  slipRate.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& slipRateSection = slipRate.section();
  assert(!slipRateSection.isNull());  

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(slipRate, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());

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
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  _slipRateVertex.resize(spaceDim);
  double_array vCoordsGlobal(spaceDim);
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    coordinates->restrictPoint(*v_iter, 
			       &vCoordsGlobal[0], vCoordsGlobal.size());
    normalizer.dimensionalize(&vCoordsGlobal[0], vCoordsGlobal.size(),
			      lengthScale);
    
    int err = _dbSlipRate->query(&_slipRateVertex[0], _slipRateVertex.size(), 
				 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbSlipRate->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_slipRateVertex[0], _slipRateVertex.size(),
				 velocityScale);

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

    slipRateSection->updatePoint(*v_iter, &_slipRateVertex[0]);
    slipTimeSection->updatePoint(*v_iter, &_slipTimeVertex);
  } // for

  // Close databases
  _dbSlipRate->close();
  _dbSlipTime->close();
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::ConstRateSlipFn::slip(topology::Field<topology::SubMesh>* slip,
				      const double t)
{ // slip
  assert(0 != slip);
  assert(0 != _parameters);

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = slip->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const topology::Field<topology::SubMesh>& slipRate = 
    _parameters->get("slip rate");
  const ALE::Obj<RealSection>& slipRateSection = slipRate.section();
  assert(!slipRateSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    slipRateSection->restrictPoint(*v_iter, &_slipRateVertex[0],
				   _slipRateVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);

    const double relTime = t - _slipTimeVertex;
    const double elapsedTime = (relTime > 0) ? relTime : 0.0;
    _slipRateVertex *= elapsedTime; // Convert slip rate to slip
    
    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipRateVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * (1 + _slipRateVertex.size()));
} // slip

// ----------------------------------------------------------------------
// Get increment of slip on fault surface between time t0 and t1.
void
pylith::faults::ConstRateSlipFn::slipIncr(
				      topology::Field<topology::SubMesh>* slip,
				      const double t0,
				      const double t1)
{ // slipIncr
  assert(0 != slip);
  assert(0 != _parameters);

  // Get vertices in fault mesh
  const ALE::Obj<SieveMesh>& sieveMesh = slip->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const label_sequence::iterator verticesBegin = vertices->begin();
  const label_sequence::iterator verticesEnd = vertices->end();

  // Get sections
  const topology::Field<topology::SubMesh>& slipRate = 
    _parameters->get("slip rate");
  const ALE::Obj<RealSection>& slipRateSection = slipRate.section();
  assert(!slipRateSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    slipRateSection->restrictPoint(*v_iter, &_slipRateVertex[0],
				   _slipRateVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);

    const double relTime0 = t0 - _slipTimeVertex;
    const double relTime1 = t1 - _slipTimeVertex;
    double elapsedTime = 0.0;
    if (relTime0 > 0)
      elapsedTime = t1 - t0;
    else if (relTime1 > 0)
      elapsedTime = t1 - _slipTimeVertex;
    _slipRateVertex *= elapsedTime; // Convert slip rate to slip
    
    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipRateVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * (4 + _slipRateVertex.size()));
} // slipIncr

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
