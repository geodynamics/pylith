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

#include "FaultCohesiveImpulses.hh" // implementation of object methods

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
pylith::faults::FaultCohesiveImpulses::FaultCohesiveImpulses(void) :
  _threshold(1.0e-6),
  _dbImpulseAmp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveImpulses::~FaultCohesiveImpulses(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveImpulses::deallocate(void)
{ // deallocate
  FaultCohesiveLagrange::deallocate();

  // :TODO: Use shared pointers for amplitudes of impulses
  _dbImpulseAmp = 0;
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for amplitudes of the impulses.
void
pylith::faults::FaultCohesiveImpulses::dbImpulseAmp(spatialdata::spatialdb::SpatialDB* db)
{ // dbImpulseAmp
  _dbImpulseAmp = db;
} // dbImpulseAmp
  
// ----------------------------------------------------------------------
// Set indices of fault degrees of freedom associated with
void
pylith::faults::FaultCohesiveImpulses::impulseDOF(const int* flags,
						  const int size)
{ // impulseDOF
  if (size > 0) {
    assert(flags);
  } // if

  _impulseDOF.resize(size);
  for (int i=0; i < size; ++i) {
    _impulseDOF[i] = flags[i];
  } // for
} // impulseDOF

// ----------------------------------------------------------------------
// Set threshold for nonzero impulse amplitude.
void
pylith::faults::FaultCohesiveImpulses::threshold(const PylithScalar value)
{ // threshold
  if (value < 0) {
    std::ostringstream msg;
    msg << "Threshold (" << value << ") for nonzero amplitudes of impulses "
      "must be nonnegative";
    throw std::runtime_error(msg.str());
  } // if

  _threshold = value;
} // threshold

// ----------------------------------------------------------------------
// Get number of impulses.
int
pylith::faults::FaultCohesiveImpulses::numImpulses(void) const
{ // numImpulses
  assert(_faultMesh);
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  MPI_Comm comm = faultSieveMesh->comm();
  int numImpulsesLocal = _impulsePoints.size();
  int numImpulses = 0;
  MPI_Allreduce(&numImpulsesLocal, &numImpulses, 1, MPI_INT, MPI_SUM, comm);

  return numImpulses;
} // numImpulses

// ----------------------------------------------------------------------
// Get number of components for impulses at each point.
int
pylith::faults::FaultCohesiveImpulses::numComponents(void) const
{ // numComponents
  return _impulseDOF.size();
} // numComponents

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveImpulses::initialize(const topology::Mesh& mesh,
						  const PylithScalar upDir[3])
{ // initialize
  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Setup impulses
  _setupImpulses();
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveImpulses::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const PylithScalar t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(fields);
  assert(_fields);
  assert(_logger);

  const int setupEvent = _logger->eventId("FaIR setup");
  _logger->eventBegin(setupEvent);

  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  dispRel.zero();
  // Set impulse corresponding to current time.
  _setRelativeDisp(dispRel, int(t+0.1));

  // Transform slip from local (fault) coordinate system to relative
  // displacement field in global coordinate system
  const topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");
  FaultCohesiveLagrange::faultToGlobal(&dispRel, orientation);

  _logger->eventEnd(setupEvent);

  FaultCohesiveLagrange::integrateResidual(residual, t, fields);

} // integrateResidual

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveImpulses::vertexField(const char* name,
						   const topology::SolutionFields* fields)
{ // vertexField
  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const topology::Field<topology::SubMesh>& orientation = _fields->get("orientation");

  PylithScalar scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(0);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("strike_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get("orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(1);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("dip_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get("orientation").section();
    assert(!orientationSection.isNull());
    const int space = (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(space);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(dirSection);
    buffer.label("normal_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strcasecmp("impulse_amplitude", name)) {
    topology::Field<topology::SubMesh>& amplitude = _fields->get("impulse amplitude");
    return amplitude;

  } else if (0 == strcasecmp("area", name)) {
    topology::Field<topology::SubMesh>& area = _fields->get("area");
    return area;

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    _calcTractionsChange(&buffer, dispT);
    return buffer;

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else


  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveImpulses::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
  return buffer;
} // vertexField

 
// ----------------------------------------------------------------------
// Setup amplitudes of impulses.
void
pylith::faults::FaultCohesiveImpulses::_setupImpulses(void)
{ // _setupImpulses
  // If no impulse amplitude specified, leave method
  if (!_dbImpulseAmp)
    return;

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  const int spaceDim = _quadrature->spaceDim();

  // Create section to hold amplitudes of impulses.
  _fields->add("impulse amplitude", "impulse_amplitude");
  topology::Field<topology::SubMesh>& amplitude = _fields->get("impulse amplitude");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  const int fiberDim = 1;
  amplitude.newSection(dispRel, fiberDim);
  amplitude.allocate();
  amplitude.scale(lengthScale);

  PylithScalar amplitudeVertex;
  const ALE::Obj<RealSection>& amplitudeSection = amplitude.section();
  assert(!amplitudeSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
  assert(cs);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  scalar_array coordsVertex(spaceDim);
  const ALE::Obj<RealSection>& coordsSection = faultSieveMesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());

  assert(_dbImpulseAmp);
  _dbImpulseAmp->open();
  const char* valueNames[1] = { "slip" };
  _dbImpulseAmp->queryVals(valueNames, 1);

  std::map<int, int> pointOrder;
  int count = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    coordsSection->restrictPoint(v_fault, &coordsVertex[0], coordsVertex.size());
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);

    amplitudeVertex = 0.0;
    int err = _dbImpulseAmp->query(&amplitudeVertex, 1, &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find amplitude for Green's function impulses at \n" << "(";
      for (int i = 0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") in fault " << label() << "\n"
	  << "using spatial database '" << _dbImpulseAmp->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if

    if (fabs(amplitudeVertex) < _threshold) {
      amplitudeVertex = 0.0;
    } // if
    _normalizer->nondimensionalize(&amplitudeVertex, 1, lengthScale);

    if (fabs(amplitudeVertex) > 0.0) {
      pointOrder[iVertex] = count;
      ++count;
    } // if

    assert(1 == amplitudeSection->getFiberDimension(v_fault));
    amplitudeSection->updatePoint(v_fault, &amplitudeVertex);
  } // for

  // Close properties database
  _dbImpulseAmp->close();

  //amplitude.view("IMPULSE AMPLITUDE"); // DEBUGGING

  _setupImpulseOrder(pointOrder);
} // _setupImpulses


// ----------------------------------------------------------------------
// Setup order of implulses.
void
pylith::faults::FaultCohesiveImpulses::_setupImpulseOrder(const std::map<int,int>& pointOrder)
{ // _setupImpulseOrder
  // Order of impulses is set by processor rank and order of points in
  // mesh, using only those points with nonzero amplitudes.

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Gather number of points on each processor.
  int numImpulsesLocal = pointOrder.size();
  const int commSize = faultSieveMesh->commSize();
  const int commRank = faultSieveMesh->commRank();
  int_array numImpulsesAll(commSize);
  MPI_Comm comm = faultSieveMesh->comm();
  PetscErrorCode err = 0;
  err = MPI_Allgather(&numImpulsesLocal, 1, MPI_INT, &numImpulsesAll[0], 1, MPI_INT, comm);CHKERRXX(err);

  int localOffset = 0;
  for (int i=0; i < commRank; ++i) {
    localOffset += numImpulsesAll[i];
  } // for

  const int ncomps = _impulseDOF.size();

  _impulsePoints.clear();
  ImpulseInfoStruct impulseInfo;
  const std::map<int,int>::const_iterator pointOrderEnd = pointOrder.end();
  for (std::map<int,int>::const_iterator piter=pointOrder.begin(); piter != pointOrderEnd; ++piter) {
    impulseInfo.indexCohesive = piter->first;
    const int offset = localOffset+piter->second;
    for (int icomp=0; icomp < ncomps; ++icomp) {
      const int impulse = ncomps*offset + icomp;
      impulseInfo.indexDOF = _impulseDOF[icomp];
      _impulsePoints[impulse] = impulseInfo;
    } // for
  } // for

#if 0 // DEBUGGING
  const ALE::Obj<RealSection>& amplitudeSection = _fields->get("impulse amplitude").section();
  assert(!amplitudeSection.isNull());
  int impulse = 0;
  for (int irank=0; irank < commSize; ++irank) {
    MPI_Barrier(comm);
    if (commRank == irank) {
      std::cout << "RANK: " << commRank << ", # impulses: " << _impulsePoints.size() << std::endl;
      const srcs_type::const_iterator impulsePointsEnd = _impulsePoints.end();
      for (srcs_type::const_iterator piter=_impulsePoints.begin(); piter != impulsePointsEnd; ++piter) {
	const int impulse = piter->first;
	const ImpulseInfoStruct& info = piter->second;
	const PylithScalar* amplitudeVertex = amplitudeSection->restrictPoint(_cohesiveVertices[info.indexCohesive].fault);
	std::cout << "["<<irank<<"]: " << impulse << " -> (" << info.indexCohesive << "," << info.indexDOF << "), v_fault: " << _cohesiveVertices[info.indexCohesive].fault << ", amplitude: " << amplitudeVertex[0] << std::endl;
      } // for
    } // if
  } // for
#endif
} // _setupImpulseOrder


// ----------------------------------------------------------------------
// Set relative displacemet associated with impulse.
void
pylith::faults::FaultCohesiveImpulses::_setRelativeDisp(const topology::Field<topology::SubMesh>& dispRel,
							const int impulse)
{ // _setRelativeDisp
  assert(_fields);

  // If no impulse amplitude specified, leave method
  if (!_dbImpulseAmp)
    return;

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  const ALE::Obj<RealSection>& amplitudeSection = _fields->get("impulse amplitude").section();
  assert(!amplitudeSection.isNull());
  
  const ALE::Obj<RealSection>& dispRelSection = dispRel.section();
  assert(!dispRelSection.isNull());

  scalar_array dispRelVertex(spaceDim);
  dispRelVertex = 0.0;
    
  const srcs_type::const_iterator& impulseInfo = _impulsePoints.find(impulse);
  if (impulseInfo != _impulsePoints.end()) {
    const int iVertex = impulseInfo->second.indexCohesive;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;

    // Get amplitude of slip impulse
    assert(1 == amplitudeSection->getFiberDimension(v_fault));
    const PylithScalar* amplitudeVertex = amplitudeSection->restrictPoint(v_fault);
    assert(amplitudeVertex);

    const int indexDOF = impulseInfo->second.indexDOF;
    assert(indexDOF >= 0 && indexDOF < spaceDim);
    dispRelVertex[indexDOF] = amplitudeVertex[0];

    assert(dispRelVertex.size() == dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updatePoint(v_fault, &dispRelVertex[0]);
  } // if

#if 0 // DEBUGGING
  std::cout << "impulse: " << impulse << std::endl;
  dispRel.view("DISP RELATIVE"); // DEBUGGING
#endif
} // _setRelativeDisp


// End of file 
