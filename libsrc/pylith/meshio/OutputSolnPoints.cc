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

#include "pylith/problems/Problem.hh" // USES Problem
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/geocoords/Converter.hh" // USES Converter
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputSolnPoints::_pyreComponent = "outputsoln";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPoints::OutputSolnPoints(pylith::problems::Problem* const problem) :
    OutputSoln(problem),
    _stationsMesh(NULL),
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

    delete _stationsMesh; _stationsMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set station names and coordinates of points .
void
pylith::meshio::OutputSolnPoints::stations(const PylithReal* stationCoords,
                                           const PylithInt numStations,
                                           const PylithInt spaceDim,
                                           const char* const* stationNames,
                                           const PylithInt numStationNames) {
    PYLITH_METHOD_BEGIN;

    assert(stationCoords && stationNames);
    assert(numStations == numStationNames);

    // Copy station coordinates.
    const PylithInt size = numStations * spaceDim;
    _stationCoords.resize(size);
    for (PylithInt i = 0; i < size; ++i) {
        _stationCoords[i] = stationCoords[i];
    } // for

    // Copy station names.
    _stationNames.resize(numStationNames);
    for (PylithInt i = 0; i < numStationNames; ++i) {
        _stationNames[i] = stationNames[i];
    } // for

    PYLITH_METHOD_END;
} // points

// ----------------------------------------------------------------------
// Setup interpolator.
void
pylith::meshio::OutputSolnPoints::_setupInterpolator(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = DMInterpolationDestroy(&_interpolator); PYLITH_CHECK_ERROR(err);
    assert(!_interpolator);

    const spatialdata::geocoords::CoordSys* csMesh = mesh.coordsys(); assert(csMesh);
    const int spaceDim = csMesh->spaceDim();

    // Setup interpolator object
    PetscDM dmMesh = mesh.dmMesh(); assert(dmMesh);

    err = DMInterpolationCreate(mesh.comm(), &_interpolator); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationSetDim(_interpolator, spaceDim); PYLITH_CHECK_ERROR(err);
    err = DMInterpolationAddPoints(_interpolator, _stationCoords.size(), (PetscReal*) &_stationCoords[0]); PYLITH_CHECK_ERROR(err);
    const PetscBool pointsAllProcs = PETSC_TRUE;
    err = DMInterpolationSetUp(_interpolator, dmMesh, pointsAllProcs); PYLITH_CHECK_ERROR(err);

    // Create mesh corresponding to local stations.
    const int meshDim = 0;
    delete _stationsMesh; _stationsMesh = new pylith::topology::Mesh(meshDim, mesh.comm()); assert(_stationsMesh);

    const int numVertices = _interpolator->n; // Number of local stations
    PylithScalar* stationCoordsLocal = NULL;
    err = VecGetArray(_interpolator->coords, &stationCoordsLocal); PYLITH_CHECK_ERROR(err);
    const int sizeLocal = numVertices * spaceDim;
    scalar_array vertices(sizeLocal);
    for (PylithInt i = 0; i < sizeLocal; ++i) {
        vertices[i] = stationCoordsLocal[i];
    } // for
    const int numCells = numVertices;
    const int numCorners = 1;
    pylith::int_array cells(numCells * numCorners);
    for (PylithInt i = 0; i < numCells; ++i) {
        cells[i] = i;
    } // for
    const bool isParallel = true;
    MeshBuilder::buildMesh(_stationsMesh, &vertices, numVertices, spaceDim,
                           cells, numCells, numCorners, meshDim, isParallel);
    err = VecRestoreArray(_interpolator->coords, &stationCoordsLocal); PYLITH_CHECK_ERROR(err);

    // Set coordinate system and create nondimensionalized coordinates
    _stationsMesh->coordsys(mesh.coordsys());

    PylithReal lengthScale = 1.0;
    err = DMPlexGetScale(mesh.dmMesh(), PETSC_UNIT_LENGTH, &lengthScale); PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(_stationsMesh->dmMesh(), PETSC_UNIT_LENGTH, lengthScale); PYLITH_CHECK_ERROR(err);

#if 0 // DEBUGGING
    _stationsMesh->view("::ascii_info_detail");
#endif

    // Upate station names to only local stations.
    pylith::string_vector stationNamesLocal(numVertices);
    const size_t numStations = _stationNames.size();
    for (PylithInt iVertex = 0; iVertex < numVertices; ++iVertex) {
        // Find point in array of all points to get index for station name.
        for (size_t iStation = 0; iStation < numStations; ++iStation) {
            const PylithReal tolerance = 1.0e-6;
            PylithReal dist = 0.0;
            for (int iDim = 0; iDim < spaceDim; ++iDim) {
                dist += pow(_stationCoords[iStation*spaceDim+iDim] - vertices[iVertex*spaceDim+iDim], 2);
            } // for
            if (sqrt(dist) < tolerance) {
                stationNamesLocal[iVertex] = _stationNames[iStation];
                break;
            } // if
        } // for
    } // for

    _stationNames = stationNamesLocal;
    _stationCoords.resize(0);

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

    PetscDM pointsDMMesh = _stationsMesh->dmMesh(); assert(pointsDMMesh);

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

    OutputSoln::appendVertexField(t, fieldInterp, *_stationsMesh);
#endif

    PYLITH_METHOD_END;
} // appendVertexField


// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::OutputSolnPoints::_writePointNames(void) {
    PYLITH_METHOD_BEGIN;

    assert(_writer);
    _writer->writePointNames(_stationNames, *_stationsMesh);

    PYLITH_METHOD_END;
} // writePointNames

// End of file
