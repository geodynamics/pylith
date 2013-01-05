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
  return _impulsePoints.size();
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
    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    return buffer;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, 0, 0, orientationVec);
    buffer.label("strike_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, 0, 1, orientationVec);
    buffer.label("dip_dir");
    buffer.scale(1.0);
    return buffer;

  } else if (0 == strcasecmp("normal_dir", name)) {
    PetscSection orientationSection = _fields->get("orientation").petscSection();
    Vec          orientationVec     = _fields->get("orientation").localVector();
    assert(orientationSection);assert(orientationVec);
    const int space = (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer = _fields->get("buffer (vector)");
    buffer.copy(orientationSection, 0, space, orientationVec);
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
  PetscErrorCode err;

  // Create section to hold amplitudes of impulses.
  _fields->add("impulse amplitude", "impulse_amplitude");
  topology::Field<topology::SubMesh>& amplitude = _fields->get("impulse amplitude");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  const int fiberDim = 1;
  amplitude.newSection(dispRel, fiberDim);
  amplitude.allocate();
  amplitude.scale(lengthScale);

  PetscSection amplitudeSection = amplitude.petscSection();
  Vec          amplitudeVec     = amplitude.localVector();
  PetscScalar *amplitudeArray;
  assert(amplitudeSection);assert(amplitudeVec);

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
  assert(cs);

  DM faultDMMesh = _faultMesh->dmMesh();
  assert(faultDMMesh);

  scalar_array coordsVertex(spaceDim);
  PetscSection coordSection;
  Vec          coordVec;
  PetscScalar *coords;
  err = DMPlexGetCoordinateSection(faultDMMesh, &coordSection);CHECK_PETSC_ERROR(err);
  err = DMGetCoordinatesLocal(faultDMMesh, &coordVec);CHECK_PETSC_ERROR(err);

  assert(_dbImpulseAmp);
  _dbImpulseAmp->open();
  const char* valueNames[1] = { "slip" };
  _dbImpulseAmp->queryVals(valueNames, 1);

  std::map<int, int> pointOrder;
  int count = 0;
  const int numVertices = _cohesiveVertices.size();
  err = VecGetArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecGetArray(amplitudeVec, &amplitudeArray);CHECK_PETSC_ERROR(err);
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    PetscInt cdof, coff;

    err = PetscSectionGetDof(coordSection, v_fault, &cdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(coordSection, v_fault, &coff);CHECK_PETSC_ERROR(err);
    assert(cdof == spaceDim);

    for(PetscInt d = 0; d < cdof; ++d) coordsVertex[d] = coords[coff+d];
    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(), lengthScale);
    PetscInt adof, aoff;

    err = PetscSectionGetDof(amplitudeSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(amplitudeSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(adof == fiberDim);

    amplitudeArray[aoff] = 0.0;
    int err = _dbImpulseAmp->query(&amplitudeArray[aoff], 1, &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find amplitude for Green's function impulses at \n" << "(";
      for (int i = 0; i < spaceDim; ++i)
        msg << "  " << coordsVertex[i];
      msg << ") in fault " << label() << "\n"
          << "using spatial database '" << _dbImpulseAmp->label() << "'.";
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
  err = VecRestoreArray(coordVec, &coords);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(amplitudeVec, &amplitudeArray);CHECK_PETSC_ERROR(err);

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

  DM faultDMMesh = _faultMesh->dmMesh();
  PetscErrorCode err;
  assert(faultDMMesh);

  // Gather number of points on each processor.
  const int numImpulsesLocal = pointOrder.size();
  MPI_Comm    comm = ((PetscObject) faultDMMesh)->comm;
  PetscMPIInt commSize, commRank;
  err = MPI_Comm_size(comm, &commSize);CHECK_PETSC_ERROR(err);
  err = MPI_Comm_rank(comm, &commRank);CHECK_PETSC_ERROR(err);
  int_array numImpulsesAll(commSize);
  err = MPI_Allgather((void *) &numImpulsesLocal, 1, MPI_INT, (void *) &numImpulsesAll[0], commSize, MPI_INT, comm);CHECK_PETSC_ERROR(err);
  
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
      for (int i=0; i < _impulsePoints.size(); ++i, ++impulse) {
	const ImpulseInfoStruct& info = _impulsePoints[impulse];
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
  PetscErrorCode err;

  PetscSection amplitudeSection = _fields->get("impulse amplitude").petscSection();
  Vec          amplitudeVec     = _fields->get("impulse amplitude").localVector();
  PetscScalar *amplitudeArray;
  assert(amplitudeSection);assert(amplitudeVec);
  
  PetscSection dispRelSection = dispRel.petscSection();
  Vec          dispRelVec     = dispRel.localVector();
  PetscScalar *dispRelArray;
  assert(dispRelSection);assert(dispRelVec);

  const srcs_type::const_iterator& impulseInfo = _impulsePoints.find(impulse);
  err = VecGetArray(amplitudeVec, &amplitudeArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);
  if (impulseInfo != _impulsePoints.end()) {
    const int iVertex = impulseInfo->second.indexCohesive;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;

    // Get amplitude of slip impulse
    PetscInt adof, aoff;

    err = PetscSectionGetDof(amplitudeSection, v_fault, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(amplitudeSection, v_fault, &aoff);CHECK_PETSC_ERROR(err);
    assert(adof == 1);
    PetscInt drdof, droff;

    err = PetscSectionGetDof(dispRelSection, v_fault, &drdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispRelSection, v_fault, &droff);CHECK_PETSC_ERROR(err);
    assert(drdof == spaceDim);
    for(PetscInt d = 0; d < drdof; ++d) {dispRelArray[droff+d] = 0.0;}
    const int indexDOF = impulseInfo->second.indexDOF;

    assert(indexDOF >= 0 && indexDOF < spaceDim);
    dispRelArray[droff+indexDOF] = amplitudeArray[aoff];
  } // if
  err = VecRestoreArray(amplitudeVec, &amplitudeArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispRelVec, &dispRelArray);CHECK_PETSC_ERROR(err);

#if 0 // DEBUGGING
  std::cout << "impulse: " << impulse << std::endl;
  dispRel.view("DISP RELATIVE"); // DEBUGGING
#endif
} // _setRelativeDisp


// End of file 
