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

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc
#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveKin::deallocate(void)
{ // deallocate
  FaultCohesiveLagrange::deallocate();

  // :TODO: Use shared pointers for earthquake sources
} // deallocate

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrcs(const char* const * names,
					 const int numNames,
					 EqKinSrc** sources,
					 const int numSources)
{ // eqsrcs
  PYLITH_METHOD_BEGIN;

  assert(numNames == numSources);

  // :TODO: Use shared pointers for earthquake sources
  _eqSrcs.clear();
  for (int i = 0; i < numSources; ++i) {
    if (0 == sources[i])
      throw std::runtime_error("Null earthquake source.");
    _eqSrcs[std::string(names[i])] = sources[i];
  } // for

  PYLITH_METHOD_END;
} // eqsrcs

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const topology::Mesh& mesh,
					     const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(src);
    src->initialize(*_faultMesh, *_normalizer);
  } // for

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateResidual(const topology::Field& residual,
						    const PylithScalar t,
						    topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("FaIR setup");
  _logger->eventBegin(setupEvent);

  topology::Field& dispRel = _fields->get("relative disp");
  dispRel.zeroAll();
  // Compute slip field at current time step
  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(src);
    if (t >= src->originTime())
      src->slip(&dispRel, t);
  } // for

  // Transform slip from local (fault) coordinate system to relative
  // displacement field in global coordinate system
  const topology::Field& orientation = _fields->get("orientation");
  FaultCohesiveLagrange::faultToGlobal(&dispRel, orientation);

  _logger->eventEnd(setupEvent);

  FaultCohesiveLagrange::integrateResidual(residual, t, fields);


  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveKin::vertexField(const char* name,
                                              const topology::SolutionFields* fields)
{ // vertexField
  PYLITH_METHOD_BEGIN;

  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);

  const int cohesiveDim = _faultMesh->dimension();

  const topology::Field& orientation = _fields->get("orientation");

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  if (0 == strcasecmp("slip", name)) {
    const topology::Field& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    PYLITH_METHOD_RETURN(buffer);

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "strike_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "dip_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("normal_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "normal_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen + 1);

    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(s_iter->second->finalSlip());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ? 
      std::string("final_slip_") + std::string(value) : "final_slip";
    buffer.label(label.c_str());

    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    const std::string value = std::string(name).substr(timeStrLen + 1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferScalarField();
    topology::Field& buffer = _fields->get("buffer (scalar)");
    buffer.copy(s_iter->second->slipTime());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ? 
      std::string("slip_time_") + std::string(value) : "slip_time";
    buffer.label(label.c_str());
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(fields);
    const topology::Field& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    _calcTractionsChange(&buffer, dispT);
    PYLITH_METHOD_RETURN(buffer);

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else


  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveKin::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field& buffer = _fields->get("buffer (vector)");
  PYLITH_METHOD_RETURN(buffer);
} // vertexField

// End of file 
