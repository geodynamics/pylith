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

#include "FaultCohesiveImpulses.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

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
  PYLITH_METHOD_BEGIN;

  FaultCohesiveLagrange::deallocate();

  _dbImpulseAmp = 0; // :TODO: Use shared pointers for amplitudes of impulses

  PYLITH_METHOD_END;
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
  if (size > 0)
    assert(flags);

  _impulseDOF.resize(size);
  for (int i=0; i < size; ++i)
    _impulseDOF[i] = flags[i];
} // impulseDOF

// ----------------------------------------------------------------------
// Set threshold for nonzero impulse amplitude.
void
pylith::faults::FaultCohesiveImpulses::threshold(const PylithScalar value)
{ // threshold
  if (value < 0) {
    std::ostringstream msg;
    msg << "Threshold (" << value << ") for nonzero amplitudes of impulses must be nonnegative";
    throw std::runtime_error(msg.str());
  } // if

  _threshold = value;
} // threshold

// ----------------------------------------------------------------------
// Get number of impulses.
int
pylith::faults::FaultCohesiveImpulses::numImpulses(void) const
{ // numImpulses
  MPI_Comm comm = _faultMesh->comm();
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
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Setup impulses
  _setupImpulses();

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term that do
// not require assembly across cells, vertices, or processors.
void
pylith::faults::FaultCohesiveImpulses::integrateResidual(const topology::Field& residual,
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
  // Set impulse corresponding to current time.
  _setRelativeDisp(dispRel, int(t+0.1));

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
pylith::faults::FaultCohesiveImpulses::vertexField(const char* name,
						   const topology::SolutionFields* fields)
{ // vertexField
  PYLITH_METHOD_BEGIN;

  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);

  const int cohesiveDim = _faultMesh->dimension();

  const topology::Field& orientation = _fields->get("orientation");

  if (0 == strcasecmp("slip", name)) {
    const topology::Field& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    buffer.complete();
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

  } else if (0 == strcasecmp("impulse_amplitude", name)) {
    topology::Field& amplitude = _fields->get("impulse amplitude");
    _allocateBufferScalarField();
    topology::Field& buffer = _fields->get("buffer (scalar)");
    buffer.copy(amplitude);
    buffer.label("impulse_amplitude");
    buffer.complete();
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("area", name)) {
    topology::Field& area = _fields->get("area");
    PYLITH_METHOD_RETURN(area);

  } else if (0 == strcasecmp("traction_change", name)) {
    assert(fields);
    const topology::Field& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    _calcTractionsChange(&buffer, dispT);
    PYLITH_METHOD_RETURN(buffer);

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '" << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else


  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveImpulses::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field& buffer = _fields->get("buffer (vector)");
  PYLITH_METHOD_RETURN(buffer);
} // vertexField

 
// ----------------------------------------------------------------------
// Setup amplitudes of impulses.
void
pylith::faults::FaultCohesiveImpulses::_setupImpulses(void)
{ // _setupImpulses
  PYLITH_METHOD_BEGIN;

  // If no impulse amplitude specified, leave method
  if (!_dbImpulseAmp)
    PYLITH_METHOD_END;

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  // Create section to hold amplitudes of impulses.
  _fields->add("impulse amplitude", "impulse_amplitude");
  topology::Field& amplitude = _fields->get("impulse amplitude");
  topology::Field& dispRel = _fields->get("relative disp");
  const int fiberDim = 1;
  amplitude.newSection(dispRel, fiberDim);
  amplitude.allocate();
  amplitude.scale(lengthScale);
  amplitude.vectorFieldType(topology::FieldBase::SCALAR);

  topology::VecVisitorMesh amplitudeVisitor(amplitude);
  PetscScalar *amplitudeArray = amplitudeVisitor.localArray();

  PetscErrorCode err;
  PetscDM amplitudeDM = amplitude.mesh().dmMesh();assert(amplitudeDM);
  PetscSection amplitudeSection = amplitude.localSection();assert(amplitudeSection);
  PetscSection amplitudeGlobalSection = NULL;
  PetscSF sf = NULL;
  err = DMGetPointSF(amplitudeDM, &sf);PYLITH_CHECK_ERROR(err);
  err = PetscSectionCreateGlobalSection(amplitudeSection, sf, PETSC_TRUE, &amplitudeGlobalSection);PYLITH_CHECK_ERROR(err);

  scalar_array coordsVertex(spaceDim);
  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);
  topology::CoordsVisitor coordsVisitor(faultDMMesh);
  PetscScalar *coordsArray = coordsVisitor.localArray();

  assert(_dbImpulseAmp);
  _dbImpulseAmp->open();
  const char* valueNames[1] = { "slip" };
  _dbImpulseAmp->queryVals(valueNames, 1);

  std::map<int, int> pointOrder;
  int count = 0;
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Skip clamped vertices
    if (v_fault < 0) {
      continue;
    } // if

    PetscInt goff;
    err = PetscSectionGetOffset(amplitudeGlobalSection, v_fault, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) {
      continue;
    } // if

    const PetscInt coff = coordsVisitor.sectionOffset(v_fault);
    assert(spaceDim == coordsVisitor.sectionDof(v_fault));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      coordsVertex[d] = coordsArray[coff+d];
    } // for
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);

    const PetscInt aoff = amplitudeVisitor.sectionOffset(v_fault);
    assert(fiberDim == amplitudeVisitor.sectionDof(v_fault));

    amplitudeArray[aoff] = 0.0;
    int err = _dbImpulseAmp->query(&amplitudeArray[aoff], 1, &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find amplitude for Green's function impulses at " << "(";
      for (int i = 0; i < spaceDim; ++i)
        msg << "  " << coordsVertex[i];
      msg << ") in fault " << label() << "using spatial database '" << _dbImpulseAmp->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if

    if (fabs(amplitudeArray[aoff]) < _threshold) {
      amplitudeArray[aoff] = 0.0;
    } // if
    _normalizer->nondimensionalize(&amplitudeArray[aoff], 1, lengthScale);

    if (fabs(amplitudeArray[aoff]) > 0.0) {
      pointOrder[iVertex] = count;
      ++count;
    } // if
  } // for
  err = PetscSectionDestroy(&amplitudeGlobalSection);PYLITH_CHECK_ERROR(err);

  // Close properties database
  _dbImpulseAmp->close();

  //amplitude.view("IMPULSE AMPLITUDE"); // DEBUGGING

  _setupImpulseOrder(pointOrder);

  PYLITH_METHOD_END;
} // _setupImpulses


// ----------------------------------------------------------------------
// Setup order of implulses.
void
pylith::faults::FaultCohesiveImpulses::_setupImpulseOrder(const std::map<int,int>& pointOrder)
{ // _setupImpulseOrder
  PYLITH_METHOD_BEGIN;

  // Order of impulses is set by processor rank and order of points in
  // mesh, using only those points with nonzero amplitudes.

  PetscDM faultDMMesh = _faultMesh->dmMesh();assert(faultDMMesh);

  // Gather number of points on each processor.
  int numImpulsesLocal = pointOrder.size();
  MPI_Comm comm = _faultMesh->comm();
  PetscMPIInt commSize, commRank;
  PetscErrorCode err;
  err = MPI_Comm_size(comm, &commSize);PYLITH_CHECK_ERROR(err);
  err = MPI_Comm_rank(comm, &commRank);PYLITH_CHECK_ERROR(err);
  int_array numImpulsesAll(commSize);
  err = MPI_Allgather(&numImpulsesLocal, 1, MPI_INT, &numImpulsesAll[0], 1, MPI_INT, comm);PYLITH_CHECK_ERROR(err);
  
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
  topology::VecVisitorMesh amplitudeVisitor(_fields->get("impulse amplitude"));
  const PetscScalar* amplitudeArray = amplitudeVisitor.localArray();
  int impulse = 0;
  for (int irank=0; irank < commSize; ++irank) {
    MPI_Barrier(comm);
    if (commRank == irank) {
      std::cout << "RANK: " << commRank << ", # impulses: " << _impulsePoints.size() << std::endl;
      const srcs_type::const_iterator impulsePointsEnd = _impulsePoints.end();
      for (srcs_type::const_iterator piter=_impulsePoints.begin(); piter != impulsePointsEnd; ++piter) {
	const int impulse = piter->first;
	const ImpulseInfoStruct& info = piter->second;
	const PetscInt aoff = amplitudeVisitor.sectionOffset(_cohesiveVertices[info.indexCohesive].fault);
	assert(1 == amplitudeVisitor.sectionDof(_cohesiveVertices[info.indexCohesive].fault));
	std::cout << "["<<irank<<"]: " << impulse << " -> (" << info.indexCohesive << "," << info.indexDOF << "), v_fault: " << _cohesiveVertices[info.indexCohesive].fault << ", amplitude: " << amplitudeArray[aoff] << std::endl;
      } // for
    } // if
  } // for
#endif


  PYLITH_METHOD_END;
} // _setupImpulseOrder


// ----------------------------------------------------------------------
// Set relative displacemet associated with impulse.
void
pylith::faults::FaultCohesiveImpulses::_setRelativeDisp(const topology::Field& dispRel,
							const int impulse)
{ // _setRelativeDisp
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  // If no impulse amplitude specified, leave method
  if (!_dbImpulseAmp)
    PYLITH_METHOD_END;

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  topology::Field& amplitude = _fields->get("impulse amplitude");
  topology::VecVisitorMesh amplitudeVisitor(amplitude);
  const PetscScalar* amplitudeArray = amplitudeVisitor.localArray();

  topology::VecVisitorMesh dispRelVisitor(dispRel);
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  const srcs_type::const_iterator& impulseInfo = _impulsePoints.find(impulse);
  if (impulseInfo != _impulsePoints.end()) {
    const int iVertex = impulseInfo->second.indexCohesive;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int e_lagrange = _cohesiveVertices[iVertex].lagrange;

    // Skip clamped vertices
    if (e_lagrange < 0) {
      PYLITH_METHOD_END;
    } // if

    // Get amplitude of slip impulse
    const PetscInt aoff = amplitudeVisitor.sectionOffset(v_fault);
    assert(1 == amplitudeVisitor.sectionDof(v_fault));

    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispRelArray[droff+d] = 0.0;
    } // for

    const int indexDOF = impulseInfo->second.indexDOF;
    assert(indexDOF >= 0 && indexDOF < spaceDim);
    dispRelArray[droff+indexDOF] = amplitudeArray[aoff];
  } // if

#if 0 // DEBUGGING
  std::cout << "impulse: " << impulse << std::endl;
  dispRel.view("DISP RELATIVE"); // DEBUGGING
#endif

  PYLITH_METHOD_END;
} // _setRelativeDisp


// End of file 
