// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/ReverseCuthillMcKee.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/error.hh" // USES PylithCallPetsc()

// ----------------------------------------------------------------------
// Reorder vertices and cells in mesh.
void
pylith::topology::ReverseCuthillMcKee::reorder(topology::Mesh* mesh) {
    assert(mesh);

    PetscDMLabel dmLabel = NULL;
    PetscDM dmOrig = mesh->getDM();
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    PylithCallPetsc(DMGetLabel(dmOrig, labelName, &dmLabel));assert(dmLabel);

    PetscIS permutation = NULL;
    PetscDM dmNew = NULL;
    PylithCallPetsc(DMPlexGetOrdering(dmOrig, MATORDERINGRCM, dmLabel, &permutation));
    PylithCallPetsc(DMPlexPermute(dmOrig, permutation, &dmNew));
    PylithCallPetsc(ISDestroy(&permutation));
    mesh->setDM(dmNew, "domain");

    // Verify that all material points (cells) are consecutive.
    PetscIS valuesIS = NULL;
    PetscInt numValues = 0;
    const PetscInt* values = NULL;
    PylithCallPetsc(DMGetLabel(dmNew, labelName, &dmLabel));assert(dmLabel);
    PylithCallPetsc(DMLabelGetValueIS(dmLabel, &valuesIS));
    PylithCallPetsc(ISGetLocalSize(valuesIS, &numValues));
    PylithCallPetsc(ISGetIndices(valuesIS, &values));
    for (PetscInt iValue = 0; iValue < numValues; ++iValue) {
        PetscIS pointsIS = NULL;
        PetscInt numPoints = 0;
        const PetscInt* points = NULL;
        PylithCallPetsc(DMLabelGetStratumIS(dmLabel, values[iValue], &pointsIS));
        PylithCallPetsc(ISGetLocalSize(pointsIS, &numPoints));
        PylithCallPetsc(ISGetIndices(pointsIS, &points));
        for (PetscInt iPoint = 1; iPoint < numPoints; ++iPoint) {
            if (points[iPoint] - points[iPoint-1] != 1) {
                // Cleanup
                PylithCallPetsc(ISRestoreIndices(pointsIS, &points));
                PylithCallPetsc(ISDestroy(&pointsIS));
                PylithCallPetsc(ISRestoreIndices(valuesIS, &values));
                PylithCallPetsc(ISDestroy(&valuesIS));

                std::ostringstream msg;
                msg << "Cells for label '" << labelName << "' with value " << values[iValue] << " are not consecutive (" << points[iPoint] << " and " << points[iPoint-1] << ").";
                throw std::runtime_error(msg.str());
            } // if
        } // for
        PylithCallPetsc(ISRestoreIndices(pointsIS, &points));
        PylithCallPetsc(ISDestroy(&pointsIS));
    } // for
    PylithCallPetsc(ISRestoreIndices(valuesIS, &values));
    PylithCallPetsc(ISDestroy(&valuesIS));
} // reorder


// End of file
