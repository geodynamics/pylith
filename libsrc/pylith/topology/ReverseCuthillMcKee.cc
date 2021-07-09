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

#include "ReverseCuthillMcKee.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Reorder vertices and cells in mesh.
void
pylith::topology::ReverseCuthillMcKee::reorder(topology::Mesh* mesh) {
    assert(mesh);
    PetscErrorCode err = 0;

    PetscDMLabel dmLabel = NULL;
    PetscDM dmOrig = mesh->getDM();
    const char* const labelName = pylith::topology::Mesh::getCellsLabelName();
    err = DMGetLabel(dmOrig, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);

    PetscIS permutation = NULL;
    PetscDM dmNew = NULL;
    err = DMPlexGetOrdering(dmOrig, MATORDERINGRCM, dmLabel, &permutation);PYLITH_CHECK_ERROR(err);
    err = DMPlexPermute(dmOrig, permutation, &dmNew);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&permutation);PYLITH_CHECK_ERROR(err);
    mesh->setDM(dmNew);

    // Verify that all material points (cells) are consecutive.
    PetscIS valuesIS = NULL;
    PetscInt numValues = 0;
    const PetscInt* values = NULL;
    err = DMGetLabel(dmNew, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetValueIS(dmLabel, &valuesIS);PYLITH_CHECK_ERROR(err);
    err = ISGetLocalSize(valuesIS, &numValues);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(valuesIS, &values);PYLITH_CHECK_ERROR(err);
    for (PetscInt iValue = 0; iValue < numValues; ++iValue) {
        PetscIS pointsIS = NULL;
        PetscInt numPoints = 0;
        const PetscInt* points = NULL;
        err = DMLabelGetStratumIS(dmLabel, values[iValue], &pointsIS);PYLITH_CHECK_ERROR(err);
        err = ISGetLocalSize(pointsIS, &numPoints);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointsIS, &points);PYLITH_CHECK_ERROR(err);
        for (PetscInt iPoint = 1; iPoint < numPoints; ++iPoint) {
            if (points[iPoint] - points[iPoint-1] != 1) {
                // Cleanup
                err = ISRestoreIndices(pointsIS, &points);PYLITH_CHECK_ERROR(err);
                err = ISDestroy(&pointsIS);PYLITH_CHECK_ERROR(err);
                err = ISRestoreIndices(valuesIS, &values);PYLITH_CHECK_ERROR(err);
                err = ISDestroy(&valuesIS);PYLITH_CHECK_ERROR(err);

                std::ostringstream msg;
                msg << "Cells for label material-id " << values[iValue] << " are not consecutive (" << points[iPoint] << " and " << points[iPoint-1] << ").";
                throw std::runtime_error(msg.str());
            } // if
        } // for
        err = ISRestoreIndices(pointsIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointsIS);PYLITH_CHECK_ERROR(err);
    } // for
    err = ISRestoreIndices(valuesIS, &values);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&valuesIS);PYLITH_CHECK_ERROR(err);
} // reorder


// End of file
