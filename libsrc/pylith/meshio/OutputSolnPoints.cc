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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputSolnPoints.hh" // implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" \
    // USES Nondimensional

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPoints::OutputSolnPoints(pylith::problems::Problem* const problem) :
    OutputSoln(problem),
    _pointsMesh(NULL),
    _interpolator(NULL)
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
pylith::meshio::OutputSolnPoints::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputSoln::deallocate();

    PetscErrorCode err = DMInterpolationDestroy(&_interpolator); PYLITH_CHECK_ERROR(err);

    delete _pointsMesh; _pointsMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set station names and coordinates of points .
void
pylith::meshio::OutputSolnPoints::points(const PylithReal* points,
                                         const int numPoints,
                                         const int spaceDim,
                                         const char* const* stationNames,
                                         const int numStations,
                                         const spatialdata::geocoords::CoordSys* coordsys) {
    // 1. Convert points to mesh coordinates system.
    //
    // 2. Nondimensionalize coordinates.
    //
    // 3. Create serial mesh.

#if 0
    assert(points);

    // Create nondimensionalized array of point locations
    const int size = numPoints * spaceDim;
    scalar_array pointsNondim(size);
    for (int i = 0; i < size; ++i) {
        pointsNondim[i] = points[i] / normalizer.lengthScale();
    } // for

    // Create serial mesh corresponding to points.
    const int meshDim = 0;
    delete _pointsMesh; _pointsMesh = new pylith::topology::Mesh(meshDim, _mesh->comm()); assert(_pointsMesh);

    const int numPointsLocal = _interpolator->n;
    PylithScalar* pointsLocal = NULL;
    err = VecGetArray(_interpolator->coords, &pointsLocal); PYLITH_CHECK_ERROR(err);
    scalar_array pointsArray(numPointsLocal*spaceDim); // Array of vertex coordinates for local mesh.
    const int sizeLocal = numPointsLocal*spaceDim;
    for (int i = 0; i < sizeLocal; ++i) {
        // Must scale by length scale because we gave interpolator nondimensioned coordinates
        pointsArray[i] = pointsLocal[i]*normalizer.lengthScale();
    } // for
    int_array cells(numPointsLocal);
    for (int i = 0; i < numPointsLocal; ++i) {
        cells[i] = i;
    } // for
    const int numCells = numPointsLocal;
    const int numCorners = 1;
    const bool isParallel = true;
    MeshBuilder::buildMesh(_pointsMesh, &pointsArray, numPointsLocal, spaceDim,
                           cells, numCells, numCorners, meshDim, isParallel);
    err = VecRestoreArray(_interpolator->coords, &pointsLocal); PYLITH_CHECK_ERROR(err);

    // Set coordinate system and create nondimensionalized coordinates
    _pointsMesh->coordsys(_mesh->coordsys());
    topology::MeshOps::nondimensionalize(_pointsMesh, normalizer);

    if (!_fields) {
        _fields = new topology::Fields(*_pointsMesh); assert(_fields);
    } // if

    // Copy station names. :TODO: Reorder to match output (pointsLocal).
    _stations.resize(numPointsLocal);
    for (int iLocal = 0; iLocal < numPointsLocal; ++iLocal) {
        // Find point in array of points to get index for station name.
        for (int iAll = 0; iAll < numPoints; ++iAll) {
            const PylithScalar tolerance = 1.0e-6;
            PylithScalar dist = 0.0;
            for (int iDim = 0; iDim < spaceDim; ++iDim) {
                dist += pow(points[iAll*spaceDim+iDim] - pointsArray[iLocal*spaceDim+iDim], 2);
            } // for
            if (sqrt(dist) < tolerance) {
                _stations[iLocal] = names[iAll];
                break;
            } // if
        } // for
    } // for
#endif
    PYLITH_METHOD_END;

} // points

// ----------------------------------------------------------------------
// Setup interpolator.
void
pylith::meshio::OutputSolnPoints::_setupInterpolator(void) {
    PYLITH_METHOD_BEGIN;

    // Make private method?

#if 0
    assert(mesh);
    assert(points);

    assert(!_interpolator); // Insure clean starting point

    _mesh = mesh;

    // Create nondimensionalized array of point locations
    const int size = numPoints * spaceDim;
    scalar_array pointsNondim(size);
    for (int i = 0; i < size; ++i) {
        pointsNondim[i] = points[i] / normalizer.lengthScale();
    } // for

#if 0 // DEBUGGING
    std::cout << "OUTPUT SOLN POINTS (dimensioned)" << std::endl;
    for (int i = 0; i < numPoints; ++i) {
        for (int iDim = 0; iDim < spaceDim; ++iDim) {
            std::cout << " " << points[i*spaceDim+iDim];
        } // for
        std::cout << "\n";
    } // for
#endif

    const spatialdata::geocoords::CoordSys* csMesh = mesh->coordsys(); assert(csMesh);
    assert(csMesh->spaceDim() == spaceDim);

    // Setup interpolator object
    PetscDM dmMesh = _mesh->dmMesh(); assert(dmMesh);
    PetscErrorCode err = 0;

    err = DMInterpolationCreate(_mesh->comm(), &_interpolator); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationSetDim(_interpolator, spaceDim); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationAddPoints(_interpolator, numPoints, (PetscReal*) &pointsNondim[0]); PYLITH_CHECK_ERROR(err);
    const PetscBool pointsAllProcs = PETSC_TRUE;
    err = DMInterpolationSetUp(_interpolator, dmMesh, pointsAllProcs); PYLITH_CHECK_ERROR(err);

    // Create mesh corresponding to points.
    const int meshDim = 0;
    delete _pointsMesh; _pointsMesh = new topology::Mesh(meshDim, _mesh->comm()); assert(_pointsMesh);

    const int numPointsLocal = _interpolator->n;
    PylithScalar* pointsLocal = NULL;
    err = VecGetArray(_interpolator->coords, &pointsLocal); PYLITH_CHECK_ERROR(err);
    scalar_array pointsArray(numPointsLocal*spaceDim); // Array of vertex coordinates for local mesh.
    const int sizeLocal = numPointsLocal*spaceDim;
    for (int i = 0; i < sizeLocal; ++i) {
        // Must scale by length scale because we gave interpolator nondimensioned coordinates
        pointsArray[i] = pointsLocal[i]*normalizer.lengthScale();
    } // for
    int_array cells(numPointsLocal);
    for (int i = 0; i < numPointsLocal; ++i) {
        cells[i] = i;
    } // for
    const int numCells = numPointsLocal;
    const int numCorners = 1;
    const bool isParallel = true;
    MeshBuilder::buildMesh(_pointsMesh, &pointsArray, numPointsLocal, spaceDim,
                           cells, numCells, numCorners, meshDim, isParallel);
    err = VecRestoreArray(_interpolator->coords, &pointsLocal); PYLITH_CHECK_ERROR(err);

    // Set coordinate system and create nondimensionalized coordinates
    _pointsMesh->coordsys(_mesh->coordsys());
    topology::MeshOps::nondimensionalize(_pointsMesh, normalizer);

#if 0 // DEBUGGING
    _pointsMesh->view("::ascii_info_detail");
#endif

    if (!_fields) {
        _fields = new topology::Fields(*_pointsMesh); assert(_fields);
    } // if

    // Copy station names. :TODO: Reorder to match output (pointsLocal).
    _stations.resize(numPointsLocal);
    for (int iLocal = 0; iLocal < numPointsLocal; ++iLocal) {
        // Find point in array of points to get index for station name.
        for (int iAll = 0; iAll < numPoints; ++iAll) {
            const PylithScalar tolerance = 1.0e-6;
            PylithScalar dist = 0.0;
            for (int iDim = 0; iDim < spaceDim; ++iDim) {
                dist += pow(points[iAll*spaceDim+iDim] - pointsArray[iLocal*spaceDim+iDim], 2);
            } // for
            if (sqrt(dist) < tolerance) {
                _stations[iLocal] = names[iAll];
                break;
            } // if
        } // for
    } // for
#endif

    PYLITH_METHOD_END;
} // setupInterpolator


// ----------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputSolnPoints::_writeDataStep(const PylithReal t,
                                                 const PylithInt tindex,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_END;
} // _writeDataSet

// ----------------------------------------------------------------------
// Append finite-element field to file.
void
pylith::meshio::OutputSolnPoints::_interpolateField(void) {
    PYLITH_METHOD_BEGIN;

#if 1
    PYLITH_COMPONENT_ERROR(":TODO: Implement OutputSolnPoints::interpolateField().");
    throw std::logic_error(":TODO: Implement OutputSolnPoints::interpolateField().");
#else
    assert(_mesh);
    assert(_fields);

    PetscDM pointsDMMesh = _pointsMesh->dmMesh(); assert(pointsDMMesh);

    PetscErrorCode err;

    PetscDM dmMesh = field.dmMesh(); assert(dmMesh);
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
    err = MPI_Allreduce(&reallocate, &reallocateGlobal, 1, MPIU_INT, MPI_LOR, fieldInterp.mesh().comm()); PYLITH_CHECK_ERROR(err);
    if (reallocateGlobal) {
        pylith::topology::Field::Discretization fe;
        fe.basisOrder = 0;
        fe.quadOrder = 0;
        fe.isBasisContinuous = true;
        fe.feSpace = pylith::topology::Field::POLYNOMIAL_SPACE;
        const pylith::string_vector& subfieldNames = field.subfieldNames();
        for (size_t i = 0; i < subfieldNames.size(); ++i) {
            const pylith::topology::Field::SubfieldInfo& sinfo = field.subfieldInfo(subfieldNames[i].c_str());
            fieldInterp.subfieldAdd(sinfo.description, fe);
        } // for
        fieldInterp.subfieldsSetup();
        fieldInterp.allocate();
    } // if

    fieldInterp.label(field.label());
    fieldInterp.zeroLocal();

    err = DMInterpolationSetDof(_interpolator, fiberDim); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationEvaluate(_interpolator, dmMesh, field.localVector(), fieldInterp.localVector()); PYLITH_CHECK_ERROR(err);

    OutputSoln::appendVertexField(t, fieldInterp, *_pointsMesh);
#endif

    PYLITH_METHOD_END;
} // appendVertexField


// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::OutputSolnPoints::_writePointNames(void) {
    PYLITH_METHOD_BEGIN;

    assert(_writer);
    _writer->writePointNames(_stations, *_pointsMesh);

    PYLITH_METHOD_END;
} // writePointNames

// End of file
