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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputSolnPoints.hh" // implementation of class methods

#include "DataWriter.hh" // USES DataWriter
#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

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
  PYLITH_METHOD_BEGIN;

  OutputManager::deallocate();

  if (_interpolator) {
    PetscErrorCode err = DMInterpolationDestroy(&_interpolator);PYLITH_CHECK_ERROR(err);
  } // if

  _mesh = 0; // :TODO: Use shared pointer
  delete _pointsMesh; _pointsMesh = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get mesh associated with points.
const pylith::topology::Mesh&
pylith::meshio::OutputSolnPoints::pointsMesh(void)
{ // pointsMesh
  PYLITH_METHOD_BEGIN;

  assert(_pointsMesh);
  PYLITH_METHOD_RETURN(*_pointsMesh);
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
  PYLITH_METHOD_BEGIN;

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

  const spatialdata::geocoords::CoordSys* csMesh = mesh->coordsys();assert(csMesh);
  assert(csMesh->spaceDim() == spaceDim);

  // Setup interpolator object
  PetscDM dmMesh = _mesh->dmMesh();assert(dmMesh);
  PetscErrorCode err = 0;

  err = DMInterpolationCreate(_mesh->comm(), &_interpolator);PYLITH_CHECK_ERROR(err);
  err = DMInterpolationSetDim(_interpolator, spaceDim);PYLITH_CHECK_ERROR(err);
  err = DMInterpolationAddPoints(_interpolator, numPoints, (PetscReal*) &pointsNondim[0]);PYLITH_CHECK_ERROR(err);
  const PetscBool pointsAllProcs = PETSC_TRUE;
  err = DMInterpolationSetUp(_interpolator, dmMesh, pointsAllProcs);PYLITH_CHECK_ERROR(err);

  // Create mesh corresponding to points.
  const int meshDim = 0;
  delete _pointsMesh; _pointsMesh = new topology::Mesh(meshDim);assert(_pointsMesh);
  topology::MeshOps::createDMMesh(_pointsMesh, meshDim, _mesh->comm(), "points");

  const int numPointsLocal = _interpolator->n;
  PylithScalar* pointsLocal = NULL;
  err = VecGetArray(_interpolator->coords, &pointsLocal);PYLITH_CHECK_ERROR(err);
  scalar_array pointsArray(numPointsLocal*spaceDim);
  const int sizeLocal = numPointsLocal*spaceDim;
  for (int i=0; i < sizeLocal; ++i) {
    // Must scale by length scale because we gave interpolator nondimensioned coordinates
    pointsArray[i] = pointsLocal[i]*normalizer.lengthScale();
  } // for
  int_array cells(numPointsLocal);
  for (int i=0; i < numPointsLocal; ++i) {
    cells[i] = i;
  } // for
  const int numCells = numPointsLocal;
  const int numCorners = 1;
  const bool interpolate = false;
  const bool isParallel = true;
  MeshBuilder::buildMesh(_pointsMesh, &pointsArray, numPointsLocal, spaceDim,
			 cells, numCells, numCorners, meshDim, interpolate, isParallel);
  err = VecRestoreArray(_interpolator->coords, &pointsLocal);PYLITH_CHECK_ERROR(err);

  // Set coordinate system and create nondimensionalized coordinates
  _pointsMesh->coordsys(_mesh->coordsys());
  topology::MeshOps::nondimensionalize(_pointsMesh, normalizer);
  
#if 0 // DEBUGGING
  _pointsMesh->view("POINTS MESH");
#endif

  if (!_fields) {
    _fields = new topology::Fields(*_pointsMesh);assert(_fields);
  } // if

  PYLITH_METHOD_END;
} // setupInterpolator


// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputSolnPoints::open(const topology::Mesh& mesh,
				       const int numTimeSteps,
				       const char* label,
				       const int labelId)
{ // open
  PYLITH_METHOD_BEGIN;

  assert(!label);
  assert(!labelId);

  assert(_pointsMesh);
  OutputManager::open(*_pointsMesh, numTimeSteps, label, labelId);

  PYLITH_METHOD_END;
} // open


// ----------------------------------------------------------------------
// Setup file for writing fields at time step.
void
pylith::meshio::OutputSolnPoints::openTimeStep(const PylithScalar t,
					       const topology::Mesh& mesh,
					       const char* label,
					       const int labelId)
{ // openTimeStep
  PYLITH_METHOD_BEGIN;

  assert(!label);
  assert(!labelId);

  assert(_pointsMesh);
  OutputManager::openTimeStep(t, *_pointsMesh, label, labelId);

  PYLITH_METHOD_END;
} // openTimeStep


// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputSolnPoints::appendVertexField(const PylithScalar t,
						    topology::Field& field,
						    const topology::Mesh& mesh)
{ // appendVertexField
  PYLITH_METHOD_BEGIN;

  assert(_mesh);
  assert(_fields);

  PetscDM pointsDMMesh = _pointsMesh->dmMesh();assert(pointsDMMesh);

  PetscErrorCode err;

  PetscDM dmMesh = field.dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt numVertices = verticesStratum.size();

  PetscInt fiberDimLocal = 0;
  if (numVertices > 0) {
    topology::VecVisitorMesh fieldVisitor(field);
    fiberDimLocal = fieldVisitor.sectionDof(vStart);
  } // if
  PetscInt fiberDim = 0;
  MPI_Allreduce(&fiberDimLocal, &fiberDim, 1, MPIU_INT, MPI_MAX, field.mesh().comm());
  assert(fiberDim > 0);

  // Create field if necessary for interpolated values or recreate
  // field if mismatch in size between buffer and field.
  const std::string& fieldName = std::string(field.label()) + " (interpolated)";
  if (!_fields->hasField(fieldName.c_str())) {
    _fields->add(fieldName.c_str(), field.label());
  } // if

  topology::Field& fieldInterp = _fields->get(fieldName.c_str());
  // The decision to reallocate a field must be collective
  PetscInt reallocate = numVertices*fiberDim != fieldInterp.sectionSize();
  PetscInt reallocateGlobal = 0;
  err = MPI_Allreduce(&reallocate, &reallocateGlobal, 1, MPIU_INT, MPI_LOR, fieldInterp.mesh().comm());PYLITH_CHECK_ERROR(err);
  if (reallocateGlobal) {
    fieldInterp.newSection(topology::FieldBase::VERTICES_FIELD, fiberDim);
    fieldInterp.allocate();
  } // if

  fieldInterp.label(field.label());
  fieldInterp.vectorFieldType(field.vectorFieldType());
  fieldInterp.scale(field.scale());
  fieldInterp.zeroAll();

  const char* context = fieldName.c_str();
  fieldInterp.createScatter(*_pointsMesh, context);

  PetscVec fieldInterpVec = fieldInterp.vector(context);assert(fieldInterpVec);
  err = DMInterpolationSetDof(_interpolator, fiberDim);PYLITH_CHECK_ERROR(err);
  err = DMInterpolationEvaluate(_interpolator, dmMesh, field.localVector(), fieldInterpVec);PYLITH_CHECK_ERROR(err);

  fieldInterp.scatterGlobalToLocal(context);

  OutputManager::appendVertexField(t, fieldInterp, *_pointsMesh);

  PYLITH_METHOD_END;
} // appendVertexField


// ----------------------------------------------------------------------
// Append finite-element cell field to file.
void
pylith::meshio::OutputSolnPoints::appendCellField(const PylithScalar t,
						  topology::Field& field,
						  const char* label,
						  const int labelId)
{ // appendCellField
  throw std::logic_error("OutputSolnPoints::appendCellField() not implemented.");
} // appendCellField

// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::OutputSolnPoints::writePointNames(const char* const* names,
						  const int numNames)
{ // writePointNames
  PYLITH_METHOD_BEGIN;

  assert(_writer);
  _writer->writePointNames(names, numNames, *_pointsMesh);

  PYLITH_METHOD_END;
} // writePointNames

// End of file 
