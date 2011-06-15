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

#include "BruneSlipFn.hh" // implementation of object methods

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
  SlipTimeFn::deallocate();

  _dbFinalSlip = 0; // :TODO: Use shared pointer.
  _dbSlipTime = 0; // :TODO: Use shared pointer.
  _dbRiseTime = 0; // :TODO: Use shared pointer.
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize slip time function.
void
pylith::faults::BruneSlipFn::initialize(
			    const topology::SubMesh& faultMesh,
			    const spatialdata::units::Nondimensional& normalizer,
			    const double originTime)
{ // initialize
  assert(0 != _dbFinalSlip);
  assert(0 != _dbSlipTime);
  assert(0 != _dbRiseTime);

  const spatialdata::geocoords::CoordSys* cs = faultMesh.coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();

  const double lengthScale = normalizer.lengthScale();
  const double timeScale = normalizer.timeScale();

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
  _parameters->add("final slip", "final_slip");
  topology::Field<topology::SubMesh>& finalSlip =
    _parameters->get("final slip");
  finalSlip.newSection(vertices, spaceDim);
  finalSlip.allocate();
  finalSlip.scale(lengthScale);
  finalSlip.vectorFieldType(topology::FieldBase::VECTOR);
  const ALE::Obj<RealSection>& finalSlipSection = finalSlip.section();
  assert(!finalSlipSection.isNull());  

  _parameters->add("slip time", "slip_time");
  topology::Field<topology::SubMesh>& slipTime = _parameters->get("slip time");
  slipTime.newSection(finalSlip, 1);
  slipTime.allocate();
  slipTime.scale(timeScale);
  slipTime.vectorFieldType(topology::FieldBase::SCALAR);
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());

  _parameters->add("rise time", "rise_time");
  topology::Field<topology::SubMesh>& riseTime = _parameters->get("rise time");
  riseTime.cloneSection(slipTime);
  riseTime.scale(timeScale);
  riseTime.vectorFieldType(topology::FieldBase::SCALAR);
  const ALE::Obj<RealSection>& riseTimeSection = riseTime.section();
  assert(!riseTimeSection.isNull());

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
      throw std::logic_error("Bad spatial dimension in BruneSlipFn.");
    } // switch

  _dbSlipTime->open();
  const char* slipTimeValues[] = {"slip-time"};
  _dbSlipTime->queryVals(slipTimeValues, 1);

  _dbRiseTime->open();
  const char* riseTimeValues[] = {"rise-time"};
  _dbRiseTime->queryVals(riseTimeValues, 1);

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
    
    int err = _dbFinalSlip->query(&_slipVertex[0], _slipVertex.size(), 
				 &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find slip rate at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbFinalSlip->label() << ".";
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

    err = _dbRiseTime->query(&_riseTimeVertex, 1, 
			     &vCoordsGlobal[0], vCoordsGlobal.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find rise time at (";
      for (int i=0; i < spaceDim; ++i)
	msg << "  " << vCoordsGlobal[i];
      msg << ") using spatial database " << _dbRiseTime->label() << ".";
      throw std::runtime_error(msg.str());
    } // if
    normalizer.nondimensionalize(&_riseTimeVertex, 1, timeScale);

    finalSlipSection->updatePoint(*v_iter, &_slipVertex[0]);
    slipTimeSection->updatePoint(*v_iter, &_slipTimeVertex);
    riseTimeSection->updatePoint(*v_iter, &_riseTimeVertex);
  } // for

  // Close databases
  _dbFinalSlip->close();
  _dbSlipTime->close();
  _dbRiseTime->close();
} // initialize

// ----------------------------------------------------------------------
// Get slip on fault surface at time t.
void
pylith::faults::BruneSlipFn::slip(topology::Field<topology::SubMesh>* slip,
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
  const topology::Field<topology::SubMesh>& finalSlip = 
    _parameters->get("final slip");
  const ALE::Obj<RealSection>& finalSlipSection = finalSlip.section();
  assert(!finalSlipSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const topology::Field<topology::SubMesh>& riseTime =
    _parameters->get("rise time");
  const ALE::Obj<RealSection>& riseTimeSection = riseTime.section();
  assert(!riseTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  const int spaceDim = _slipVertex.size();
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    finalSlipSection->restrictPoint(*v_iter, &_slipVertex[0],
				   _slipVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);
    riseTimeSection->restrictPoint(*v_iter, &_riseTimeVertex, 1);

    double finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += _slipVertex[i]*_slipVertex[i];
    finalSlipMag = sqrt(finalSlipMag);

    const double slip = _slipFn(t-_slipTimeVertex, finalSlipMag,
				_riseTimeVertex);
    const double scale = finalSlipMag > 0.0 ? slip / finalSlipMag : 0.0;
    _slipVertex *= scale;
    
    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * (2+8 + 3*_slipVertex.size()));
} // slip

// ----------------------------------------------------------------------
// Get increment of slip on fault surface between time t0 and t1.
void
pylith::faults::BruneSlipFn::slipIncr(
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
  const topology::Field<topology::SubMesh>& finalSlip = 
    _parameters->get("final slip");
  const ALE::Obj<RealSection>& finalSlipSection = finalSlip.section();
  assert(!finalSlipSection.isNull());
  const topology::Field<topology::SubMesh>& slipTime =
    _parameters->get("slip time");
  const ALE::Obj<RealSection>& slipTimeSection = slipTime.section();
  assert(!slipTimeSection.isNull());
  const topology::Field<topology::SubMesh>& riseTime =
    _parameters->get("rise time");
  const ALE::Obj<RealSection>& riseTimeSection = riseTime.section();
  assert(!riseTimeSection.isNull());
  const ALE::Obj<RealSection>& slipSection = slip->section();
  assert(!slipSection.isNull());

  const int spaceDim = _slipVertex.size();
  for (label_sequence::iterator v_iter=verticesBegin;
       v_iter != verticesEnd;
       ++v_iter) {
    finalSlipSection->restrictPoint(*v_iter, &_slipVertex[0],
				   _slipVertex.size());
    slipTimeSection->restrictPoint(*v_iter, &_slipTimeVertex, 1);
    riseTimeSection->restrictPoint(*v_iter, &_riseTimeVertex, 1);

    double finalSlipMag = 0.0;
    for (int i=0; i < spaceDim; ++i)
      finalSlipMag += _slipVertex[i]*_slipVertex[i];
    finalSlipMag = sqrt(finalSlipMag);

    const double slip0 = _slipFn(t0-_slipTimeVertex, finalSlipMag,
				 _riseTimeVertex);
    const double slip1 = _slipFn(t1-_slipTimeVertex, finalSlipMag,
				 _riseTimeVertex);
    const double scale = finalSlipMag > 0.0 ? 
      (slip1 - slip0) / finalSlipMag : 0.0;
    _slipVertex *= scale;

    
    // Update field
    slipSection->updateAddPoint(*v_iter, &_slipVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * (3+2*8 + 3*_slipVertex.size()));
} // slipIncr

// ----------------------------------------------------------------------
// Get final slip.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::BruneSlipFn::finalSlip(void)
{ // finalSlip
  return _parameters->get("final slip");
} // finalSlip

// ----------------------------------------------------------------------
// Get time when slip begins at each point.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::BruneSlipFn::slipTime(void)
{ // slipTime
  return _parameters->get("slip time");
} // slipTime


// End of file 
