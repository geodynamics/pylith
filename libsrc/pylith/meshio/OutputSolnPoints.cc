// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputSolnPoints.hh" // implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/geocoords/Converter.hh" // USES Converter
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPoints::OutputSolnPoints(void) :
    _pointMesh(NULL),
    _pointSoln(NULL),
    _interpolator(NULL) {
    PyreComponent::setName("outputsolnpoints");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPoints::~OutputSolnPoints(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnPoints::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputSoln::deallocate();

    if (_interpolator) {
        PetscErrorCode err = DMInterpolationDestroy(&_interpolator);PYLITH_CHECK_ERROR(err);
    } // if

    delete _pointMesh;_pointMesh = NULL;
    delete _pointSoln;_pointSoln = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set point names and coordinates of points .
void
pylith::meshio::OutputSolnPoints::setPoints(const PylithReal* pointCoords,
                                            const int numPoints,
                                            const int spaceDim,
                                            const char* const* pointNames,
                                            const int numPointNames) {
    PYLITH_METHOD_BEGIN;

    assert(pointCoords && pointNames);
    assert(numPoints == numPointNames);

    // Copy point coordinates.
    const PylithInt size = numPoints * spaceDim;
    _pointCoords.resize(size);
    for (PylithInt i = 0; i < size; ++i) {
        _pointCoords[i] = pointCoords[i];
    } // for

    // Copy point names.
    _pointNames.resize(numPointNames);
    for (PylithInt i = 0; i < numPointNames; ++i) {
        _pointNames[i] = pointNames[i];
    } // for

    PYLITH_METHOD_END;
} // setPoints


// ------------------------------------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputSolnPoints::_writeSolnStep(const PylithReal t,
                                                 const PylithInt tindex,
                                                 const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_writeSolnStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.getLabel()<<")");

    if (!_interpolator) {
        _setupInterpolator(solution);
    } // if
    assert(_pointMesh);
    assert(_pointSoln);
    _interpolateField(solution);

    const bool writePointNames = !_writer->isOpen();
    _openSolnStep(t, *_pointMesh);
    if (writePointNames) { _writePointNames(); }

    const pylith::string_vector& subfieldNames = pylith::topology::FieldOps::getSubfieldNamesDomain(solution);

    const size_t numSubfieldNames = subfieldNames.size();
    for (size_t iField = 0; iField < numSubfieldNames; iField++) {
        OutputSubfield* subfield = NULL;
        subfield = this->_getSubfield(*_pointSoln, *_pointMesh, subfieldNames[iField].c_str());assert(subfield);

        const pylith::topology::Field::SubfieldInfo& info = solution.getSubfieldInfo(subfieldNames[iField].c_str());
        subfield->extractSubfield(*_pointSoln, info.index);

        OutputObserver::_appendField(t, *subfield);
    } // for
    _closeSolnStep();

    PYLITH_METHOD_END;
} // _writeSolnStep


// ------------------------------------------------------------------------------------------------
// Get output subfield, creating if necessary.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputSolnPoints::_getSubfield(const pylith::topology::Field& field,
                                               const pylith::topology::Mesh& submesh,
                                               const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_getSubfield(field="<<field.getLabel()<<", name="<<name<<", submesh="<<typeid(submesh).name()<<")");

    if (_subfields.count(name) == 0) {
        _subfields[name] = OutputSubfield::create(field, submesh, name);
    } // if

    PYLITH_METHOD_RETURN(_subfields[name]);
}


// ------------------------------------------------------------------------------------------------
// Setup interpolator.
void
pylith::meshio::OutputSolnPoints::_setupInterpolator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    if (_interpolator) {
        err = DMInterpolationDestroy(&_interpolator);PYLITH_CHECK_ERROR(err);
    } // if
    assert(!_interpolator);

    const spatialdata::geocoords::CoordSys* csMesh = solution.getMesh().getCoordSys();assert(csMesh);
    const int spaceDim = csMesh->getSpaceDim();

    MPI_Comm comm = solution.getMesh().getComm();

    // Setup interpolator object
    PetscDM dmSoln = solution.getDM();assert(dmSoln);

    err = DMInterpolationCreate(comm, &_interpolator);PYLITH_CHECK_ERROR(err);
    err = DMInterpolationSetDim(_interpolator, spaceDim);PYLITH_CHECK_ERROR(err);
    err = DMInterpolationAddPoints(_interpolator, _pointNames.size(), (PetscReal*) &_pointCoords[0]);PYLITH_CHECK_ERROR(err);
    const PetscBool pointsAllProcs = PETSC_TRUE;
    const PetscBool ignoreOutsideDomain = PETSC_FALSE;
    err = DMInterpolationSetUp(_interpolator, dmSoln, pointsAllProcs, ignoreOutsideDomain);PYLITH_CHECK_ERROR(err);

    // Create mesh corresponding to local points.
    const size_t numPointsLocal = _interpolator->n;
    PylithScalar* pointsLocal = NULL;
    err = VecGetArray(_interpolator->coords, &pointsLocal);PYLITH_CHECK_ERROR(err);

    PylithReal lengthScale = 1.0;
    err = DMPlexGetScale(dmSoln, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);

    const spatialdata::geocoords::CoordSys* cs = solution.getMesh().getCoordSys();
    delete _pointMesh;_pointMesh = pylith::topology::MeshOps::createFromPoints(
        pointsLocal, numPointsLocal, cs, lengthScale, comm);

    // Upate point names to only local points.
    pylith::string_vector pointNamesLocal(numPointsLocal);
    const size_t numPoints = _pointNames.size();
    for (size_t iPointLocal = 0; iPointLocal < numPointsLocal; ++iPointLocal) {
        // Find point in array of all points to get index for point name.
        for (size_t iPoint = 0; iPoint < numPoints; ++iPoint) {
            const PylithReal tolerance = 1.0e-6;
            PylithReal dist = 0.0;
            for (int iDim = 0; iDim < spaceDim; ++iDim) {
                dist += pow(_pointCoords[iPoint*spaceDim+iDim] - pointsLocal[iPointLocal*spaceDim+iDim], 2);
            } // for
            if (sqrt(dist) < tolerance) {
                pointNamesLocal[iPointLocal] = _pointNames[iPoint];
                break;
            } // if
        } // for
    } // for
    err = VecGetArray(_interpolator->coords, &pointsLocal);PYLITH_CHECK_ERROR(err);

    _pointNames = pointNamesLocal;
    _pointCoords.resize(0);

    // Determine size of interpolated field that we will have.
    PetscInt numDof = 0;
    const pylith::string_vector& subfieldNames = solution.getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    for (size_t i = 0; i < numSubfields; ++i) {
        numDof += solution.getSubfieldInfo(subfieldNames[i].c_str()).description.numComponents;
    } // for
    err = DMInterpolationSetDof(_interpolator, numDof);PYLITH_CHECK_ERROR(err);

    delete _pointSoln;_pointSoln = new pylith::topology::Field(*_pointMesh);
    for (size_t i = 0; i < subfieldNames.size(); ++i) {
        const pylith::topology::Field::SubfieldInfo& sinfo = solution.getSubfieldInfo(subfieldNames[i].c_str());
        pylith::topology::Field::Discretization discretization = sinfo.fe;
        discretization.dimension = spaceDim;
        _pointSoln->subfieldAdd(sinfo.description, discretization);
    } // for
    _pointSoln->subfieldsSetup();
    _pointSoln->createDiscretization();
    _pointSoln->setLabel(solution.getLabel());
    _pointSoln->allocate();

    PYLITH_METHOD_END;
} // setupInterpolator


// ------------------------------------------------------------------------------------------------
// Append finite-element field to file.
void
pylith::meshio::OutputSolnPoints::_interpolateField(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    assert(_pointSoln);

    PetscErrorCode err;
    err = DMInterpolationEvaluate(_interpolator, solution.getDM(), solution.getLocalVector(),
                                  _pointSoln->getLocalVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // appendVertexField


// ------------------------------------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::OutputSolnPoints::_writePointNames(void) {
    PYLITH_METHOD_BEGIN;

    assert(_writer);
    _writer->writePointNames(_pointNames, *_pointMesh);

    PYLITH_METHOD_END;
} // writePointNames


// End of file
