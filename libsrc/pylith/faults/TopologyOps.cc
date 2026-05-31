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

#include "pylith/faults/TopologyOps.hh" // implementation of object methods

#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/MeshOps.hh" // USES isCohesiveCell()
#include "pylith/utils/error.hh" // USES PylithCallPetsc()
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <iostream> // USES std::cout
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::updateCohesiveLabel(const pylith::topology::Mesh* mesh,
                                                 const char* labelName,
                                                 const int labelValue) {
    assert(mesh);

    const PetscDM dmMesh = mesh->getDM();
    PetscDMLabel dmLabel = nullptr;
    PylithCallPetsc(DMGetLabel(dmMesh, labelName, &dmLabel));

    // Set label value of points to dimension (vertices=0, edges=1, faces=2)
    PetscIS pointIS = nullptr;
    const PylithInt* points = nullptr;
    PylithInt numPoints = 0;
    PylithCallPetsc(DMLabelGetStratumIS(dmLabel, labelValue, &pointIS));
    if (pointIS) {
        PylithCallPetsc(ISGetIndices(pointIS, &points));
        PylithCallPetsc(ISGetSize(pointIS, &numPoints));
    } // if
    for (PylithInt iPoint = 0; iPoint < numPoints; ++iPoint) {
        const PylithInt point = points[iPoint];

        // Clear existing point in label (need label values to be depth).
        PylithCallPetsc(DMLabelClearValue(dmLabel, point, labelValue));

        PylithInt *closure = nullptr;
        PylithInt closureSize = 0, pointDepth = 0;

        PylithCallPetsc(DMPlexGetPointDepth(dmMesh, point, &pointDepth));
        PylithCallPetsc(DMLabelSetValue(dmLabel, point, pointDepth));
        PylithCallPetsc(DMPlexGetTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure));
        for (PylithInt cl = 0; cl < closureSize * 2; cl += 2) {
            PylithCallPetsc(DMPlexGetPointDepth(dmMesh, closure[cl], &pointDepth));
            PylithCallPetsc(DMLabelSetValue(dmLabel, closure[cl], pointDepth));
        } // for
        PylithCallPetsc(DMPlexRestoreTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure));
    } // for
    if (pointIS) {
        PylithCallPetsc(ISRestoreIndices(pointIS, &points));
        PylithCallPetsc(ISDestroy(&pointIS));
    } // if

    PylithCallPetsc(DMLabelDestroyIndex(dmLabel)); // :KLUDGE).
    PylithCallPetsc(DMPlexOrientLabel(dmMesh, dmLabel));
    {
        IS valueIS;
        const PetscInt *values;
        PetscInt depth, Nv;

        PylithCallPetsc(DMPlexGetDepth(dmMesh, &depth));
        PylithCallPetsc(DMLabelGetValueIS(dmLabel, &valueIS));
        PylithCallPetsc(ISGetLocalSize(valueIS, &Nv));
        PylithCallPetsc(ISGetIndices(valueIS, &values));
        PylithCallPetsc(DMPlexRebalanceSharedLabelPoints(dmMesh, dmLabel, Nv, values, depth - 1));
        PylithCallPetsc(ISRestoreIndices(valueIS, &values));
        PylithCallPetsc(ISDestroy(&valueIS));
    }
    PylithCallPetsc(DMPlexCheckOrientationLabel(dmMesh, dmLabel));
    PylithCallPetsc(DMPlexLabelCohesiveComplete(dmMesh, dmLabel, nullptr, 0, PETSC_FALSE, NULL));
    DMLabel bdlabel = nullptr;
    PylithCallPetsc(DMGetLabel(dmMesh, "fault_edges", &bdlabel));
    PylithCallPetsc(DMPlexLabelCohesiveCheck(dmMesh, dmLabel, bdlabel));
} // updateCohesiveLabel


// ------------------------------------------------------------------------------------------------
// Set label and label value for newly created cohesive cells.
void
pylith::faults::TopologyOps::setCohesiveCellLabel(PetscDM dmMeshNew,
                                                  PetscDM dmMesh,
                                                  const char* labelName,
                                                  const PylithInt labelValue) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Stratum oldStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PylithInt pStart = oldStratum.end();
    pylith::topology::Stratum newStratum(dmMeshNew, pylith::topology::Stratum::HEIGHT, 0);
    const PylithInt pEnd = newStratum.end();
    const PylithInt cohesiveLabelValue = labelValue;

    PetscDMLabel dmLabel = NULL;
    PylithCallPetsc(DMGetLabel(dmMeshNew, labelName, &dmLabel));

    PetscSF sf = NULL;
    const PetscInt *leaves = NULL;
    PetscInt numLeaves = 0;
    PylithCallPetsc(DMGetPointSF(dmMeshNew, &sf));
    PylithCallPetsc(PetscSFGetGraph(sf, NULL, &numLeaves, &leaves, NULL));
    if (leaves) {
        for (PylithInt point = pStart; point < pEnd; ++point) {
            PetscInt loc;
            PylithCallPetsc(PetscFindInt(point, numLeaves, leaves, &loc));
            if (loc < 0) { // not a ghost cell
                DMLabelSetValue(dmLabel, point, cohesiveLabelValue);
            } // if
        } // for
    } else {
        for (PylithInt point = pStart; point < pEnd; ++point) {
            DMLabelSetValue(dmLabel, point, cohesiveLabelValue);
        } // for
    } // if/else

    PYLITH_METHOD_END;
} // setCohesiveCellLabel


// ------------------------------------------------------------------------------------------------
// Remove cohesive points from labels.
void
pylith::faults::TopologyOps::labelsRemoveCohesivePoints(PetscDM dmMeshNew) {
    PYLITH_METHOD_BEGIN;

    PylithInt numLabels = 0;
    PylithCallPetsc(DMGetNumLabels(dmMeshNew, &numLabels));

    const std::string& materialLabelName = pylith::topology::Mesh::cells_label_name;

    for (int iLabel = 0; iLabel < numLabels; ++iLabel) {
        const char* labelStr = NULL;
        PylithCallPetsc(DMGetLabelName(dmMeshNew, iLabel, &labelStr));
        const std::string labelName = std::string(labelStr);

        if ((labelName != std::string("depth"))
            && (labelName != std::string("celltype"))
            && (labelName != materialLabelName)) {
            PetscDMLabel dmLabel = PETSC_NULLPTR;
            PylithInt pStart = -1, pEnd = -1;
            PylithCallPetsc(DMGetLabel(dmMeshNew, labelStr, &dmLabel));
            PylithCallPetsc(DMLabelGetBounds(dmLabel, &pStart, &pEnd));

            PetscIS valuesIS = PETSC_NULLPTR;
            PylithCallPetsc(DMLabelGetNonEmptyStratumValuesIS(dmLabel, &valuesIS));
            PylithInt numValues = 0;
            PylithCallPetsc(ISGetLocalSize(valuesIS, &numValues));
            const PylithInt* valuesIndices = PETSC_NULLPTR;
            PylithCallPetsc(ISGetIndices(valuesIS, &valuesIndices));
            for (PylithInt point = pStart; point < pEnd; ++point) {
                DMPolytopeType ct;
                PylithCallPetsc(DMPlexGetCellType(dmMeshNew, point, &ct));
                if ((ct == DM_POLYTOPE_POINT_PRISM_TENSOR) ||
                    (ct == DM_POLYTOPE_SEG_PRISM_TENSOR) ||
                    (ct == DM_POLYTOPE_TRI_PRISM_TENSOR) ||
                    (ct == DM_POLYTOPE_QUAD_PRISM_TENSOR)) {
                    for (PylithInt iValue = 0; iValue < numValues; ++iValue) {
                        PylithCallPetsc(DMLabelClearValue(dmLabel, point, valuesIndices[iValue]));
                    } // for
                } // if
            } // for
            PylithCallPetsc(ISRestoreIndices(valuesIS, &valuesIndices));
            PylithCallPetsc(ISDestroy(&valuesIS));
        } // if
    } // for

    PYLITH_METHOD_END;
} // labelsRemoveCohesivePoints


// ------------------------------------------------------------------------------------------------
// Create buried edge label.
void
pylith::faults::TopologyOps::createBuriedEdgeLabel(PetscDM dmMeshNew,
                                                   PetscDM dmMesh,
                                                   const char* buriedEdgeLabelName,
                                                   const PylithInt buriedEdgeLabelValue,
                                                   PetscDMLabel surfaceLabel,
                                                   DMPlexTransform transform) {
    PYLITH_METHOD_BEGIN;

    bool hasBuriedEdge = false;
    PetscDMLabel buriedEdgeLabel = NULL;
    PylithInt dimMesh = 0;

    PylithCallPetsc(DMCreateLabel(dmMeshNew, buriedEdgeLabelName));
    PylithCallPetsc(DMGetLabel(dmMeshNew, buriedEdgeLabelName, &buriedEdgeLabel));
    PylithCallPetsc(DMGetDimension(dmMesh, &dimMesh));

    PylithInt pStart, pEnd;
    PylithCallPetsc(DMPlexGetChart(dmMesh, &pStart, &pEnd));
    for (PylithInt point = pStart; point < pEnd; ++point) {
        DMPolytopeType cellType;
        PetscInt value;
        PylithCallPetsc(DMPlexGetCellType(dmMesh, point, &cellType));
        PylithCallPetsc(DMLabelGetValue(surfaceLabel, point, &value));
        if (value >= 200) { // buried edge (no splitting)
            PylithInt numCellTypes; // Number of cell types.
            DMPolytopeType* newCellTypes = NULL; // List of new cell types [numCellTypes]
            PylithInt* newCellTypesSize = NULL; // Sizes for newCellTypes [numCellTypes];
            PylithInt* newPointsCones = NULL; // Cone of every point it makes.
            PylithInt* coneOrientations = NULL; // Orientation of every cone point
            PylithInt dim, pointEdge = 0;

            hasBuriedEdge = true;
            dim = DMPolytopeTypeGetDim(cellType);
            PylithCallPetsc(DMPlexTransformCellTransform(transform, cellType, point, NULL, &numCellTypes, &newCellTypes, &newCellTypesSize, &newPointsCones, &coneOrientations));
            for (PylithInt iCellType = 0; iCellType < numCellTypes; ++iCellType) {
                if ((DMPolytopeTypeGetDim(newCellTypes[iCellType]) != dim) || (dim != dimMesh-2)) {
                    continue;
                } // if
                const PetscInt cellTypeSize = newCellTypesSize[iCellType];
                assert(1 == cellTypeSize);
                const PetscInt iPoint = 0; // single point in cell type size
                PylithCallPetsc(DMPlexTransformGetTargetPoint(transform, cellType, newCellTypes[iCellType], point, iPoint, &pointEdge));
                PylithCallPetsc(DMLabelSetValue(buriedEdgeLabel, pointEdge, 1));
            } // for
        } // if

    } // for
    const int hasBuriedEdgeValueLocal = hasBuriedEdge ? 1 : 0;
    int hasBuriedEdgeValue = 0;
    int mpierr = MPI_Allreduce(&hasBuriedEdgeValueLocal, &hasBuriedEdgeValue, 1, MPI_INT, MPI_MAX, PetscObjectComm((PetscObject) dmMeshNew));assert(MPI_SUCCESS == mpierr);
    if (!hasBuriedEdgeValue) {
        PylithCallPetsc(DMRemoveLabel(dmMeshNew, buriedEdgeLabelName, NULL));
    } // if

    PYLITH_METHOD_END;
} // createBuriedEdgeLabel


// ------------------------------------------------------------------------------------------------
// Create fault mesh from cohesive cells.
void
pylith::faults::TopologyOps::createFaultFromCohesiveCells(pylith::topology::Mesh* faultMesh,
                                                          const pylith::topology::Mesh& mesh,
                                                          const char* cohesiveLabelName,
                                                          const int cohesiveLabelValue,
                                                          const char* surfaceLabel) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);

    faultMesh->setCoordSys(mesh.getCoordSys());

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscDM dmFaultMesh = NULL;

    const char* negativeLabelName = "fault_cohesive_negative_sides";
    PetscDMLabel negativeLabel = nullptr;
    const PylithInt negativeLabelValue = 1;
    { // Create label over negative sides of cohesive cells
        PylithCallPetsc(DMLabelCreate(mesh.getComm(), negativeLabelName, &negativeLabel));

        PetscDMLabel cohesiveLabel = nullptr;
        PylithCallPetsc(DMGetLabel(dmDomain, cohesiveLabelName, &cohesiveLabel));

        PetscIS pointIS = nullptr;
        const PylithInt* points = nullptr;
        PylithInt numPoints = 0;
        PylithCallPetsc(DMLabelGetStratumIS(cohesiveLabel, cohesiveLabelValue, &pointIS));
        if (pointIS) {
            PylithCallPetsc(ISGetIndices(pointIS, &points));
            PylithCallPetsc(ISGetSize(pointIS, &numPoints));
        } // if
        for (PylithInt iPoint = 0; iPoint < numPoints; ++iPoint) {
            const PylithInt point = points[iPoint];
            const PylithInt *cone = nullptr;
            PylithInt coneSize = 0;
            PylithCallPetsc(DMPlexGetConeSize(dmDomain, point, &coneSize));assert(coneSize > 0);
            PylithCallPetsc(DMPlexGetCone(dmDomain, point, &cone));
            const PylithInt negativeFace = cone[0];
            PylithCallPetsc(DMLabelSetValue(negativeLabel, negativeFace, negativeLabelValue));
        } // for
        if (pointIS) {PylithCallPetsc(ISRestoreIndices(pointIS, &points));}
        PylithCallPetsc(ISDestroy(&pointIS));
        PylithCallPetsc(DMPlexLabelComplete(dmDomain, negativeLabel));
    } // Create label over negative sides of cohesive cells

    const PetscBool markedFaces = PETSC_TRUE;
    PylithCallPetsc(DMPlexCreateSubmesh(dmDomain, negativeLabel, negativeLabelValue, markedFaces, &dmFaultMesh));
    PylithCallPetsc(DMLabelDestroy(&negativeLabel));
    PylithCallPetsc(DMPlexOrient(dmFaultMesh)); // :TODO: Is this necessary?

    PetscReal lengthScale = 1.0;
    PylithCallPetsc(DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(DMPlexSetScale(dmFaultMesh, PETSC_UNIT_LENGTH, lengthScale));

    faultMesh->setDM(dmFaultMesh, surfaceLabel);
    pylith::topology::MeshOps::checkTopology(*faultMesh);

    PYLITH_METHOD_END;
} // createFaultFromCohesiveCells


// ------------------------------------------------------------------------------------------------
// Get name of PETSc DM label for interfaces.
const char*
pylith::faults::TopologyOps::getInterfacesLabelName(void) {
    return "cohesive interface";
} // getInterfacesLabelName


// ------------------------------------------------------------------------------------------------
// Get PETSc DM label for interfaces, creating if necessary.
PetscDMLabel
pylith::faults::TopologyOps::getInterfacesLabel(PetscDM dm) {
    PYLITH_METHOD_BEGIN;
    PetscDMLabel interfacesLabel = NULL;

    const char* interfacesLabelName = TopologyOps::getInterfacesLabelName();
    PetscBool hasInterfacesLabel = PETSC_FALSE;
    if (DMHasLabel(dm, interfacesLabelName, &hasInterfacesLabel)) {
        PylithCallPetsc(DMGetLabel(dm, interfacesLabelName, &interfacesLabel));
    } else {
        PylithInt dim = 0;
        PylithInt pStart = 0;
        PylithInt pEnd = 0;
        PylithInt pMax = 0;

        PylithCallPetsc(DMGetDimension(dm, &dim));
        PylithCallPetsc(DMCreateLabel(dm, interfacesLabelName));
        PylithCallPetsc(DMGetLabel(dm, interfacesLabelName, &interfacesLabel));
        for (PylithInt iDim = 0; iDim <= dim; ++iDim) {
            PylithCallPetsc(DMPlexGetHeightStratum(dm, iDim, &pStart, &pEnd));
            PylithCallPetsc(DMPlexGetSimplexOrBoxCells(dm, iDim, NULL, &pMax));
            for (PylithInt p = pMax; p < pEnd; ++p) {
                PylithCallPetsc(DMLabelSetValue(interfacesLabel, p, 1));
            } // for
        } // for
    } // else
    assert(interfacesLabel);

    PYLITH_METHOD_RETURN(interfacesLabel);
} // getInterfacesLabel


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::getAdjacentCells(PylithInt* adjacentCellNegative,
                                              PylithInt* adjacentCellPositive,
                                              PetscDM dmDomain,
                                              const PylithInt cohesiveCell) {
    PYLITH_METHOD_BEGIN;

    const PylithInt* cone = NULL;
    PylithCallPetsc(DMPlexGetCone(dmDomain, cohesiveCell, &cone));
    PylithInt adjacentCells[2];
    for (PylithInt iCone = 0; iCone < 2; ++iCone) {
        const PylithInt* support = NULL;
        PylithInt supportSize = 0;

        PylithCallPetsc(DMPlexGetSupport(dmDomain, cone[iCone], &support));
        PylithCallPetsc(DMPlexGetSupportSize(dmDomain, cone[iCone], &supportSize));
        if (2 != supportSize) {
            PYLITH_ERROR(pylith::TopologyError, pylith::journal::logic,
                         "Inconsistent topology. Expected support of size 2 for face "
                         << cone[iCone] << " of cohesive cell " << cohesiveCell
                         <<". Support has size "<<supportSize<<".");
        } // if
        assert(2 == supportSize);
        if ((cohesiveCell != support[0]) && (cohesiveCell != support[1]) ) {
            PYLITH_ERROR(pylith::TopologyError, pylith::journal::logic,
                         "Inconsistent topology. Cohesive cell "
                         <<cohesiveCell<<" not in support of its own cone. "
                         <<"Support: "<<support[0]<< ", "<<support[1]<<".");
        } // if
        adjacentCells[iCone] = (support[0] == cohesiveCell) ? support[1] : support[0];
    } // for

    if (adjacentCellNegative) { *adjacentCellNegative = adjacentCells[0]; }
    if (adjacentCellPositive) { *adjacentCellPositive = adjacentCells[1]; }
    PYLITH_METHOD_END;
} // getAdjacentCells


// End of file
