// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "Source.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::Source::Source(void) :
    _description(""),
    _labelName(""),
    _labelValue(1) {
    //
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::Source::~Source(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::Source::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::Physics::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set descriptive label of source.
void
pylith::sources::Source::setDescription(const char* value) {
    PYLITH_COMPONENT_DEBUG("setDescription(value="<<value<<")");

    _description = value;
} // setDescription


// ------------------------------------------------------------------------------------------------
// Get label of source.
const char*
pylith::sources::Source::getDescription(void) const {
    return _description.c_str();
} // getDescription


// ------------------------------------------------------------------------------------------------
// Set name of solution subfield associated with boundary condition.
void
pylith::sources::Source::setSubfieldName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSubfieldName(value="<<value<<")");

    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        msg << "Empty string given for name of solution subfield for source '"
            << getLabelName() <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ------------------------------------------------------------------------------------------------
// Get name of solution subfield associated with boundary condition.
const char*
pylith::sources::Source::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::sources::Source::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<") empty method");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Set point names and coordinates of points .
void
pylith::sources::Source::setPoints(const PylithReal* pointCoords,
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
// Convert cartesian positions to a labeled source
void
pylith::sources::Source::locateSource(const pylith::topology::Field& solution) {
    printf("In MomentTensorForce begin\n");
    DMView(solution.getDM(), NULL);
    PetscErrorCode err;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    // transform points of source to mesh coordinates in python
    // DM from solution
    Vec vecPoints;
    DMLabel label;
    PetscSF sfPoints = NULL;
    const PetscInt *localPoints;
    const PetscSFNode *remotePoints;
    PetscInt numRoots = -1, numLeaves, dim, d;
    PetscMPIInt rank;
    PetscScalar       *a;

    err = DMGetCoordinateDim(dmSoln, &dim);PYLITH_CHECK_ERROR(err);
    err = VecCreateSeqWithArray(PETSC_COMM_SELF, dim, _pointCoords.size(),
                                &_pointCoords[0], &vecPoints);PYLITH_CHECK_ERROR(err);

    // Debug
    // PetscPrintf(PetscObjectComm((PetscObject) dmSoln), "_pointCoords\n");
    // PetscPrintf(PetscObjectComm((PetscObject) dmSoln), " x = %g\n", _pointCoords[0]);
    // PetscPrintf(PetscObjectComm((PetscObject) dmSoln), " y = %g\n", _pointCoords[1]);
    // Erzatz from ex17
    // err = VecCreateSeq(PETSC_COMM_SELF, dim, &vecPoints);PYLITH_CHECK_ERROR(err);
    // err = VecSetBlockSize(vecPoints, _pointCoords.size());PYLITH_CHECK_ERROR(err);
    // err = VecGetArray(vecPoints, &a);PYLITH_CHECK_ERROR(err);
    // for (d = 0; d < _pointCoords.size(); ++d) {
    //     a[d] = _pointCoords[d];
    // }
    // err = VecRestoreArray(vecPoints, &a);PYLITH_CHECK_ERROR(err);

    err = DMLocatePoints(dmSoln, vecPoints, DM_POINTLOCATION_NONE, &sfPoints);PYLITH_CHECK_ERROR(err);
    err = PetscSFView(sfPoints, NULL);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&vecPoints);PYLITH_CHECK_ERROR(err);
    err = DMCreateLabel(dmSoln,getLabelName());PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmSoln,getLabelName(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscSFGetGraph(sfPoints, &numRoots, &numLeaves, &localPoints, &remotePoints);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmSoln), &rank);PYLITH_CHECK_ERROR(err);
    // Debug
    // PetscPrintf(PetscObjectComm((PetscObject) dmSoln), "localPoints: %D\n", numLeaves);
    for (PetscInt p = 0; p < numLeaves; ++p) {
        if (remotePoints[p].rank == rank) {
            err = DMLabelSetValue(label, remotePoints[p].index, 2);PYLITH_CHECK_ERROR(err);
        }
    } // for
    err = PetscSFDestroy(&sfPoints);PYLITH_CHECK_ERROR(err);
    printf("In MomentTensorForce end\n");
    DMView(dmSoln, NULL);
}


// End of file
