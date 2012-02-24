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

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc
#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

//#define PRECOMPUTE_GEOMETRY
//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

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
  assert(numNames == numSources);

  // :TODO: Use shared pointers for earthquake sources
  _eqSrcs.clear();
  for (int i = 0; i < numSources; ++i) {
    if (0 == sources[i])
      throw std::runtime_error("Null earthquake source.");
    _eqSrcs[std::string(names[i])] = sources[i];
  } // for
} // eqsrcs

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const topology::Mesh& mesh,
					     const double upDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    src->initialize(*_faultMesh, *_normalizer);
  } // for
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  const int setupEvent = _logger->eventId("FaIR setup");
  _logger->eventBegin(setupEvent);

  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  dispRel.zero();
  // Compute slip field at current time step
  const srcs_type::const_iterator srcsEnd = _eqSrcs.end();
  for (srcs_type::iterator s_iter = _eqSrcs.begin(); s_iter != srcsEnd; ++s_iter) {
    EqKinSrc* src = s_iter->second;
    assert(0 != src);
    if (t >= src->originTime())
      src->slip(&dispRel, t);
  } // for

  // Transform slip from local (fault) coordinate system to relative
  // displacement field in global coordinate system
  _faultToGlobal(&dispRel);

  _logger->eventEnd(setupEvent);

  FaultCohesiveLagrange::integrateResidual(residual, t, fields);

} // integrateResidual

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveKin::vertexField(const char* name,
                                              const topology::SolutionFields* fields)
{ // vertexField
  assert(0 != _faultMesh);
  assert(0 != _quadrature);
  assert(0 != _normalizer);
  assert(0 != _fields);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const int slipStrLen = strlen("final_slip");
  const int timeStrLen = strlen("slip_time");

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    _globalToFault(&buffer);
    return buffer;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      0);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("strike_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      1);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("dip_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const int space = (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      space);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("normal_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strncasecmp("final_slip_X", name, slipStrLen)) {
    const std::string value = std::string(name).substr(slipStrLen + 1);

    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(s_iter->second->finalSlip());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ? 
      std::string("final_slip_") + std::string(value) : "final_slip";
    buffer.label(label.c_str());

    return buffer;

  } else if (0 == strncasecmp("slip_time_X", name, timeStrLen)) {
    const std::string value = std::string(name).substr(timeStrLen + 1);
    const srcs_type::const_iterator s_iter = _eqSrcs.find(value);
    assert(s_iter != _eqSrcs.end());

    // Need to append name of rupture to final slip label. Because
    // Field is const, we use a buffer.
    _allocateBufferScalarField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (scalar)");
    buffer.copy(s_iter->second->slipTime());
    assert(value.length() > 0);
    const std::string& label = (_eqSrcs.size() > 1) ? 
      std::string("slip_time_") + std::string(value) : "slip_time";
    buffer.label(label.c_str());

    return buffer;

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractionsChange(&buffer, dispT);
    return buffer;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else


  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveKin::vertexField().");

  // Satisfy return values
  assert(0 != _fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get(
    "buffer (vector)");
  return buffer;
} // vertexField

// End of file 
