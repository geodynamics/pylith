// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputSolnPoints.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "MeshBuilder.hh" // USES MeshBuilder

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPoints::OutputSolnPoints(void) :
  _mesh(0),
  _pointsMesh(0),
  _interpolator(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPoints::~OutputSolnPoints(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnPoints::deallocate(void)
{ // deallocate
  OutputManager<topology::Mesh, topology::Field<topology::Mesh> >::deallocate();

  if (_interpolator) {
    PetscErrorCode err = err = DMInterpolationDestroy(&_interpolator);CHECK_PETSC_ERROR(err);
  } // if

  _mesh = 0; // :TODO: Use shared pointer
  delete _pointsMesh; _pointsMesh = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get mesh associated with points.
const pylith::topology::Mesh&
pylith::meshio::OutputSolnPoints::pointsMesh(void)
{ // pointsMesh
  assert(_pointsMesh);
  return *_pointsMesh;
} // pointsMesh


// ----------------------------------------------------------------------
// Setup interpolator.
void
pylith::meshio::OutputSolnPoints::setupInterpolator(topology::Mesh* mesh,
						    const PylithScalar* points,
						    const int numPoints,
						    const int spaceDim,
						    const spatialdata::units::Nondimensional& normalizer)
{ // setupInterpolator
  assert(mesh);
  assert(points);

  assert(!_interpolator); // Insure clean starting point

  _mesh = mesh;

  // Create nondimensionalized array of point locations
  const int size = numPoints * spaceDim;
  scalar_array pointsNondim(size);
  for (int i=0; i < size; ++i) {
    pointsNondim[i] = points[i] / normalizer.lengthScale();
  } // for

#if 0 // DEBUGGING
  std::cout << "OUTPUT SOLN POINTS (dimensioned)" << std::endl;
  for (int i=0; i < numPoints; ++i) {
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      std::cout << " " << points[i*spaceDim+iDim];
    } // for
    std::cout << "\n";
  } // for
#endif

  const spatialdata::geocoords::CoordSys* csMesh = mesh->coordsys();
  assert(csMesh);
  assert(csMesh->spaceDim() == spaceDim);

  // Setup interpolator object
  DM dm = _mesh->dmMesh();
  MPI_Comm       comm;
  PetscErrorCode err = 0;

  assert(dm);
  err = PetscObjectGetComm((PetscObject) dm, &comm);CHECK_PETSC_ERROR(err);
  err = DMInterpolationCreate(comm, &_interpolator);CHECK_PETSC_ERROR(err);
  
  err = DMInterpolationSetDim(_interpolator, spaceDim);CHECK_PETSC_ERROR(err);

  err = DMInterpolationAddPoints(_interpolator, numPoints, (PetscReal*) &pointsNondim[0]);CHECK_PETSC_ERROR(err);
  const PetscBool pointsAllProcs = PETSC_TRUE;
  err = DMInterpolationSetUp(_interpolator, dm, pointsAllProcs);CHECK_PETSC_ERROR(err);

  // Create mesh corresponding to points.
  const int meshDim = 0;
  delete _pointsMesh; _pointsMesh = new topology::Mesh(meshDim);
  _pointsMesh->createDMMesh(meshDim);
  assert(_pointsMesh);

  const int numPointsLocal = _interpolator->n;
  PylithScalar* pointsLocal = 0;
  err = VecGetArray(_interpolator->coords, &pointsLocal);CHECK_PETSC_ERROR(err);
  scalar_array pointsArray(numPointsLocal*spaceDim);
  const int sizeLocal = numPointsLocal*spaceDim;
  for (int i=0; i < sizeLocal; ++i) {
    // Must scale by length scale because we gave interpolator nondimensioned coordinates
    pointsArray[i] = pointsLocal[i]*normalizer.lengthScale();
  } // for
  int_array cells(numPointsLocal);
  for (int i=0; i < numPointsLocal; ++i)
    cells[i] = i;
  const int numCells = numPointsLocal;
  const int numCorners = 1;
  const bool interpolate = false;
  const bool isParallel = true;
  MeshBuilder::buildMesh(_pointsMesh, &pointsArray, numPointsLocal, spaceDim,
			 cells, numCells, numCorners, meshDim, interpolate, isParallel);
  err = VecRestoreArray(_interpolator->coords, &pointsLocal);CHECK_PETSC_ERROR(err);

  // Set coordinate system and create nondimensionalized coordinates
  _pointsMesh->coordsys(_mesh->coordsys());
  _pointsMesh->nondimensionalize(normalizer);

#if 0 // DEBUGGING
  _pointsMesh->view("POINTS MESH");
#endif

  if (!_fields) {
    _fields = new topology::Fields<topology::Field<topology::Mesh> >(*_pointsMesh);
  } // if

} // setupInterpolator


// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputSolnPoints::open(const topology::Mesh& mesh,
				       const int numTimeSteps,
				       const char* label,
				       const int labelId)
{ // open
  assert(!label);
  assert(!labelId);

  assert(_pointsMesh);
  OutputManager<topology::Mesh, topology::Field<topology::Mesh> >::open(*_pointsMesh, numTimeSteps, label, labelId);
} // open


// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputSolnPoints::openTimeStep(const PylithScalar t,
					       const topology::Mesh& mesh,
					       const char* label,
					       const int labelId)
{ // openTimeStep
  assert(!label);
  assert(!labelId);

  assert(_pointsMesh);
  OutputManager<topology::Mesh, topology::Field<topology::Mesh> >::openTimeStep(t, *_pointsMesh, label, labelId);
} // openTimeStep


// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputSolnPoints::appendVertexField(const PylithScalar t,
			       topology::Field<topology::Mesh>& field,
			       const topology::Mesh& mesh)
{ // appendVertexField
  assert(_mesh);
  assert(_fields);

  DM             pointsDMMesh = _pointsMesh->dmMesh();
  PetscInt       pvStart, pvEnd;
  PetscErrorCode err;

  assert(pointsDMMesh);
  err = DMPlexGetDepthStratum(pointsDMMesh, 0, &pvStart, &pvEnd);CHECK_PETSC_ERROR(err);

  DM       dmMesh = field.dmMesh();
  PetscInt vStart, vEnd;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  PetscSection section = field.petscSection();
  PetscInt fiberDimLocal = 0;
  PetscInt fiberDim = 0;

  if (vEnd > vStart) {err = PetscSectionGetDof(section, vStart, &fiberDimLocal);CHECK_PETSC_ERROR(err);}
  MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPIU_INT, MPI_MAX, field.mesh().comm());
  assert(fiberDim > 0);

  // Create field if necessary for interpolated values or recreate
  // field if mismatch in size between buffer and field.
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Output");
  ostringstream fieldName;
  const char   *context = field.label();
  fieldName << context << " (interpolated)" << std::endl;

  if (_fields->hasField(fieldName.str().c_str())) {
    ostringstream msg;
    msg << "Field " << fieldName << "already present in manager" << std::endl;
    throw std::logic_error(msg.str());
  } // if
  _fields->add(fieldName.str().c_str(), field.label());

  topology::Field<topology::Mesh>& fieldInterp = _fields->get(fieldName.str().c_str());
  fieldInterp.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
  fieldInterp.allocate();
  fieldInterp.zero();
  fieldInterp.label(field.label());
  fieldInterp.vectorFieldType(field.vectorFieldType());
  fieldInterp.scale(field.scale());
  logger.stagePop();

  fieldInterp.createScatter(*_pointsMesh, context);
  PetscVec fieldInterpVec = fieldInterp.vector(context);
  assert(fieldInterpVec);
  
  err = DMInterpolationSetDof(_interpolator, fiberDim);CHECK_PETSC_ERROR(err);
  err = DMInterpolationEvaluate(_interpolator, dmMesh, field.localVector(), fieldInterpVec);CHECK_PETSC_ERROR(err);

  fieldInterp.scatterVectorToSection(context);

  OutputManager<topology::Mesh, topology::Field<topology::Mesh> >::appendVertexField(t, fieldInterp, *_pointsMesh);
} // appendVertexField


// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputSolnPoints::appendCellField(const PylithScalar t,
			   topology::Field<topology::Mesh>& field,
			   const char* label,
			   const int labelId)
{ // appendCellField
  throw std::logic_error("OutputSolnPoints::appendCellField() not implemented.");
} // appendCellField


// End of file 
