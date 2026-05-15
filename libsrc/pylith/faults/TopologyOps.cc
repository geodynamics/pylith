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

#include <iostream> // USES std::cout
#include <cassert> // USES assert()

#if 0
// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::createFault(pylith::topology::Mesh* faultMesh,
                                         const pylith::topology::Mesh& mesh,
                                         PetscDMLabel surfaceLabel,
                                         const int surfaceLabelValue) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);

    faultMesh->setCoordSys(mesh.getCoordSys());
    PetscDM dmDomain = mesh.getDM();assert(dmDomain);

    // Convert fault to a DM
    PetscDM dmFault = PETSC_NULLPTR;
    const char *groupName = "";

    if (surfaceLabel) {
        PylithCallPetsc(PetscObjectGetName((PetscObject)surfaceLabel, &groupName));
    } // if

    PylithInt labelHasVertices = 0;
    { // TEMPORARY: Continue to support creating lower dimension meshes using labels with vertices.
        PetscIS labelIS = NULL;
        PylithInt labelHasVerticesLocal = 0;
        if (surfaceLabel) {
            const PetscInt* labelPoints = NULL;
            PetscInt numPoints = 0;
            PylithCallPetsc(DMGetStratumIS(dmDomain, groupName, surfaceLabelValue, &labelIS));
            if (labelIS) {
                PylithCallPetsc(ISGetIndices(labelIS, &labelPoints));
                PylithCallPetsc(DMGetStratumSize(dmDomain, groupName, surfaceLabelValue, &numPoints));
            } // if

            pylith::topology::Stratum verticesStratum(dmDomain, pylith::topology::Stratum::DEPTH, 0);
            PylithInt vStart = 0, vEnd = 0;
            vStart = verticesStratum.begin();
            vEnd = verticesStratum.end();
            for (PylithInt iPoint = 0; iPoint < numPoints; ++iPoint) {
                if ((labelPoints[iPoint] >= vStart) && (labelPoints[iPoint] < vEnd) ) {
                    labelHasVerticesLocal = 1;
                    break;
                } // if
            } // if
            if (labelIS) {
                PylithCallPetsc(ISRestoreIndices(labelIS, &labelPoints));
                PylithCallPetsc(ISDestroy(&labelIS));
            } // if
        } // if
        PylithCallPetsc(MPI_Allreduce(&labelHasVerticesLocal, &labelHasVertices, 1, MPI_INT, MPI_MAX,
                                      PetscObjectComm((PetscObject) dmDomain)));

        if (labelHasVertices) {
            pythia::journal::warning_t warning("deprecated");
            warning << pythia::journal::at(__HERE__)
                    << "DEPRECATION: Creating fault mesh from label with vertices. "
                    << "This feature will be removed in v6.0. "
                    << "In the future, you will need to mark boundaries not vertices for boundary conditions."
                    << pythia::journal::endl; \
        } // if
    } // TEMPORARY

    PetscDMLabel surfaceLabelFull = PETSC_NULLPTR;
    if (surfaceLabel) {
        PylithCallPetsc(DMLabelDuplicate(surfaceLabel, &surfaceLabelFull));
        PylithCallPetsc(DMPlexLabelComplete(dmDomain, surfaceLabelFull));
    } // if

    const PetscBool markedFaces = !labelHasVertices ? PETSC_TRUE : PETSC_FALSE;
    PylithCallPetsc(DMPlexCreateSubmesh(dmDomain, surfaceLabelFull, surfaceLabelValue, markedFaces, &dmFault));
    PylithCallPetsc(DMLabelDestroy(&surfaceLabelFull));

    PylithInt maxConeSizeLocal = 0, maxConeSize = 0;
    PylithCallPetsc(DMPlexGetMaxSizes(dmFault, &maxConeSizeLocal, NULL));
    PylithCallPetsc(MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmFault)));

    if (maxConeSize <= 0) {
        PylithCallPetsc(DMDestroy(&dmFault));
        std::ostringstream msg;
        msg << "Error while creating fault. Fault '" << groupName << "' with label value "
            << surfaceLabelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Check that no cells have all vertices on the fault
    if (surfaceLabel) {
        PetscIS subpointIS;
        const PylithInt *dmpoints;
        PylithInt defaultValue, cStart, cEnd, vStart, vEnd;

        PylithCallPetsc(DMLabelGetDefaultValue(surfaceLabel, &defaultValue));
        PylithCallPetsc(DMPlexGetSubpointIS(dmFault, &subpointIS));
        PylithCallPetsc(DMPlexGetHeightStratum(dmFault, 0, &cStart, &cEnd));
        PylithCallPetsc(DMPlexGetDepthStratum(dmDomain, 0, &vStart, &vEnd));
        PylithCallPetsc(ISGetIndices(subpointIS, &dmpoints));
        for (PylithInt c = cStart; c < cEnd; ++c) {
            PetscBool invalidCell = PETSC_TRUE;
            PylithInt *closure = NULL;
            PylithInt closureSize;

            PylithCallPetsc(DMPlexGetTransitiveClosure(dmDomain, dmpoints[c], PETSC_TRUE, &closureSize, &closure));
            for (PylithInt cl = 0; cl < closureSize*2; cl += 2) {
                PylithInt value = 0;

                if ((closure[cl] < vStart) || (closure[cl] >= vEnd)) { continue;}
                PylithCallPetsc(DMLabelGetValue(surfaceLabel, closure[cl], &value));
                if (value == defaultValue) {invalidCell = PETSC_FALSE;break;}
            } // for
            PylithCallPetsc(DMPlexRestoreTransitiveClosure(dmDomain, dmpoints[c], PETSC_TRUE, &closureSize, &closure));
            if (invalidCell) {
                std::ostringstream msg;
                msg << "Ambiguous fault surface. Cell "<<dmpoints[c]<<" has all of its vertices on the fault.";
                PylithCallPetsc(ISRestoreIndices(subpointIS, &dmpoints));
                PylithCallPetsc(DMDestroy(&dmFault));
                throw std::runtime_error(msg.str());
            } // if
        } // for
        PylithCallPetsc(ISRestoreIndices(subpointIS, &dmpoints));
    } // if
    PylithCallPetsc(DMPlexOrient(dmFault));

    std::string submeshLabel = "fault_" + std::string(groupName);
    faultMesh->setDM(dmFault, submeshLabel.c_str());

    PYLITH_METHOD_END;
} // createFault


// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::create(pylith::topology::Mesh* mesh,
                                    const pylith::topology::Mesh& faultMesh,
                                    PetscDMLabel faultBdLabel,
                                    const int faultBdLabelValue,
                                    const int cohesiveLabelValue) {
    assert(mesh);
    PetscDM sdm = NULL;
    PetscDM dm = mesh->getDM();assert(dm);
    PetscDMLabel subpointMap = NULL, label = NULL, mlabel = NULL;
    PylithInt dim, cMax, cStart, cEnd, numCohesiveCellsOld;

    // Have to remember the old number of cohesive cells
    PylithCallPetsc(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    cMax = cStart;
    for (PylithInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
        if (pylith::topology::MeshOps::isCohesiveCell(dm, cell)) { break; }
    } // for
    numCohesiveCellsOld = cEnd - cMax;
    // Create cohesive cells
    PylithCallPetsc(DMPlexGetSubpointMap(faultMesh.getDM(), &subpointMap));
    PylithCallPetsc(DMLabelDuplicate(subpointMap, &label));
    PylithCallPetsc(DMLabelClearStratum(label, mesh->getDimension()));
    // Fix over-aggressive completion of boundary label
    PylithCallPetsc(DMGetDimension(dm, &dim));
    if (faultBdLabel && (dim > 2)) {
        PetscIS bdIS;
        const PylithInt *bd;
        PylithInt fStart, fEnd, n, i;

        PylithCallPetsc(DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd));
        PylithCallPetsc(DMLabelGetStratumIS(faultBdLabel, faultBdLabelValue, &bdIS));
        PylithCallPetsc(ISGetLocalSize(bdIS, &n));
        PylithCallPetsc(ISGetIndices(bdIS, &bd));
        for (i = 0; i < n; ++i) {
            const PylithInt p = bd[i];

            // Remove faces
            if ((p >= fStart) && (p < fEnd)) {
                const PylithInt *edges,   *verts, *supportA, *supportB;
                PylithInt numEdges, numVerts, supportSizeA, sA, supportSizeB, sB, val, bval, e, s;
                PetscBool found = PETSC_FALSE;

                PylithCallPetsc(DMLabelClearValue(faultBdLabel, p, faultBdLabelValue));
                // Remove the cross edge
                PylithCallPetsc(DMPlexGetCone(dm, p, &edges));
                PylithCallPetsc(DMPlexGetConeSize(dm, p, &numEdges));
                if (numEdges != 3) {
                    std::ostringstream msg;
                    msg << "Internal error while creating fault mesh. Face "<<p<<" has "<<numEdges<<" edges != 3.";
                    throw std::logic_error(msg.str());
                }
                for (e = 0; e < numEdges; ++e) {
                    PylithCallPetsc(DMPlexGetCone(dm, edges[e], &verts));
                    PylithCallPetsc(DMPlexGetConeSize(dm, edges[e], &numVerts));
                    if (numVerts != 2) {
                        std::ostringstream msg;
                        msg << "Internal error while creating fault mesh. Edge "<<edges[e]<<" has "<<numVerts<<" vertices != 2.";
                        throw std::logic_error(msg.str());
                    }
                    PylithCallPetsc(DMPlexGetSupportSize(dm, verts[0], &supportSizeA));
                    PylithCallPetsc(DMPlexGetSupport(dm, verts[0], &supportA));
                    for (s = 0, sA = 0; s < supportSizeA; ++s) {
                        PylithCallPetsc(DMLabelGetValue(label, supportA[s], &val));
                        PylithCallPetsc(DMLabelGetValue(faultBdLabel, supportA[s], &bval));
                        if (( val >= 0) && ( bval >= 0) ) { ++sA;}
                    }
                    PylithCallPetsc(DMPlexGetSupportSize(dm, verts[1], &supportSizeB));
                    PylithCallPetsc(DMPlexGetSupport(dm, verts[1], &supportB));
                    for (s = 0, sB = 0; s < supportSizeB; ++s) {
                        PylithCallPetsc(DMLabelGetValue(label, supportB[s], &val));
                        PylithCallPetsc(DMLabelGetValue(faultBdLabel, supportB[s], &bval));
                        if (( val >= 0) && ( bval >= 0) ) { ++sB;}
                    }
                    if ((sA > 2) && (sB > 2)) {
                        PylithCallPetsc(DMLabelClearValue(faultBdLabel, edges[e], faultBdLabelValue));
                        found = PETSC_TRUE;
                        break;
                    }
                }
                if (!found) {
                    std::ostringstream msg;
                    msg << "Internal error while creating fault mesh. Face "<<p<<" has no cross edge.";
                    throw std::logic_error(msg.str());
                }
            }
        }
        PylithCallPetsc(ISRestoreIndices(bdIS, &bd));
        PylithCallPetsc(ISDestroy(&bdIS));
    }
    // Completes the set of cells scheduled to be replaced
    PylithCallPetsc(DMPlexOrientLabel(dm, label));
    PylithCallPetsc(DMPlexLabelCohesiveComplete(dm, label, faultBdLabel, faultBdLabelValue, PETSC_FALSE, PETSC_FALSE, faultMesh.getDM()));
    PylithCallPetsc(DMPlexConstructCohesiveCells(dm, label, NULL, &sdm));

    const char* interfaceLabelName = pylith::topology::Mesh::cells_label_name;
    PylithCallPetsc(DMGetLabel(sdm, interfaceLabelName, &mlabel));
    if (mlabel) {
        PylithCallPetsc(DMPlexGetHeightStratum(sdm, 0, &cStart, &cEnd));
        cMax = cStart;
        for (PylithInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
            if (pylith::topology::MeshOps::isCohesiveCell(sdm, cell)) { break; }
        }
        assert(cStart == cEnd || cEnd > cMax + numCohesiveCellsOld);
        for (PylithInt cell = cMax; cell < cEnd - numCohesiveCellsOld; ++cell) {
            PylithInt onBd;

            /* Eliminate hybrid cells on the boundary of the split from cohesive label,
             * they are marked with -(cell number) since the hybrid cell number aliases vertices in the old mesh */
            PylithCallPetsc(DMLabelGetValue(label, -cell, &onBd));
            // if (onBd == dim) continue;
            PylithCallPetsc(DMLabelSetValue(mlabel, cell, cohesiveLabelValue));
        }
    }
    PylithCallPetsc(DMLabelDestroy(&label));

    PetscReal lengthScale = 1.0;
    PylithCallPetsc(DMPlexGetScale(dm, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(DMPlexSetScale(sdm, PETSC_UNIT_LENGTH, lengthScale));
    PylithCallPetsc(DMViewFromOptions(sdm, NULL, "-pylith_cohesive_dm_view"));
    mesh->setDM(sdm, "domain");
} // create


// ------------------------------------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::TopologyOps::createFaultParallel(pylith::topology::Mesh* faultMesh,
                                                 const pylith::topology::Mesh& mesh,
                                                 const int labelValue,
                                                 const char* labelName,
                                                 const char* surfaceLabel) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);

    faultMesh->setCoordSys(mesh.getCoordSys());

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscDM dmFaultMesh = NULL;

    const PetscBool hasLagrangeConstraints = PETSC_TRUE;
    PylithCallPetsc(DMPlexCreateCohesiveSubmesh(dmDomain, hasLagrangeConstraints, labelName, labelValue, &dmFaultMesh));
    PylithCallPetsc(DMViewFromOptions(dmFaultMesh, NULL, "-pylith_fault_dm_view"));
    PylithCallPetsc(DMPlexOrient(dmFaultMesh));
    std::string meshLabel = std::string("fault_") + std::string(surfaceLabel);

    PetscReal lengthScale = 1.0;
    PylithCallPetsc(DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(DMPlexSetScale(dmFaultMesh, PETSC_UNIT_LENGTH, lengthScale));

    faultMesh->setDM(dmFaultMesh, meshLabel.c_str());

    PYLITH_METHOD_END;
} // createFaultParallel


// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::classifyCellsDM(PetscDM dmDomain,
                                             PylithInt vertex,
                                             const int depth,
                                             const int faceSize,
                                             PylithInt firstCohesiveCell,
                                             PointSet& replaceCells,
                                             PointSet& noReplaceCells,
                                             const int debug) {
    // Replace all cells on a given side of the fault with a vertex on the fault
    PointSet vReplaceCells;
    PointSet vNoReplaceCells;
    const PylithInt *support;
    PylithInt supportSize, s, classifyTotal = 0;
    PetscBool modified = PETSC_FALSE;

    if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
    PylithCallPetsc(DMPlexGetSupportSize(dmDomain, vertex, &supportSize));
    PylithCallPetsc(DMPlexGetSupport(dmDomain, vertex, &support));
    for (s = 0; s < supportSize; ++s) {
        const PylithInt point = support[s];

        if (point >= firstCohesiveCell) { return;}
        if (replaceCells.find(point)   != replaceCells.end()) {vReplaceCells.insert(point);}
        if (noReplaceCells.find(point) != noReplaceCells.end()) { vNoReplaceCells.insert(point);}
        modified = PETSC_TRUE;
        ++classifyTotal;
    }
    PylithInt classifySize = vReplaceCells.size() + vNoReplaceCells.size();

    while (modified && (classifySize < classifyTotal)) {
        modified = PETSC_FALSE;
        for (s = 0; s < supportSize; ++s) {
            const PylithInt point = support[s];
            PetscBool classified = PETSC_FALSE;

            if (debug) {
                const PylithInt *cone;
                PylithInt coneSize;

                std::cout << "Checking neighbor " << vertex << std::endl;
                PylithCallPetsc(DMPlexGetConeSize(dmDomain, vertex, &coneSize));
                PylithCallPetsc(DMPlexGetCone(dmDomain, vertex, &cone));
                for (PetscInt c = 0; c < coneSize; ++c) {
                    std::cout << "  cone point " << cone[c] << std::endl;
                }
            }
            if (vReplaceCells.find(point) != vReplaceCells.end()) {
                if (debug) { std::cout << "  already in replaceCells" << std::endl;}
                continue;
            } // if
            if (vNoReplaceCells.find(point) != vNoReplaceCells.end()) {
                if (debug) { std::cout << "  already in noReplaceCells" << std::endl;}
                continue;
            } // if
            if (point >= firstCohesiveCell) {
                if (debug) { std::cout << "  already a cohesive cell" << std::endl;}
                continue;
            } // if
              // If neighbor shares a face with anyone in replaceCells, then add
            for (PointSet::const_iterator c_iter = vReplaceCells.begin(); c_iter != vReplaceCells.end(); ++c_iter) {
                const PylithInt *coveringPoints;
                PylithInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                PylithCallPetsc(DMPlexGetMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints));
                PylithCallPetsc(DMPlexRestoreMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints));
                if (numCoveringPoints == faceSize) {
                    if (debug) { std::cout << "    Scheduling " << point << " for replacement" << std::endl;}
                    vReplaceCells.insert(point);
                    modified = PETSC_TRUE;
                    classified = PETSC_TRUE;
                    break;
                } // if
            } // for
            if (classified) { continue;}
            // It is unclear whether taking out the noReplace cells will speed this up
            for (PointSet::const_iterator c_iter = vNoReplaceCells.begin(); c_iter != vNoReplaceCells.end(); ++c_iter) {
                const PylithInt *coveringPoints;
                PylithInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                PylithCallPetsc(DMPlexGetMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints));
                PylithCallPetsc(DMPlexRestoreMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints));
                if (numCoveringPoints == faceSize) {
                    if (debug) { std::cout << "    Scheduling " << point << " for no replacement" << std::endl;}
                    vNoReplaceCells.insert(point);
                    modified = PETSC_TRUE;
                    classified = PETSC_TRUE;
                    break;
                } // for
            } // for
        }
        if (debug) {
            std::cout << "classifySize: " << classifySize << std::endl;
            std::cout << "classifyTotal: " << classifyTotal << std::endl;
            std::cout << "vReplaceCells.size: " << vReplaceCells.size() << std::endl;
            std::cout << "vNoReplaceCells.size: " << vNoReplaceCells.size() << std::endl;
        }
        assert(size_t(classifySize) < vReplaceCells.size() + vNoReplaceCells.size());
        classifySize = vReplaceCells.size() + vNoReplaceCells.size();
        if (classifySize > classifyTotal) {
            std::ostringstream msg;
            msg << "Internal error classifying cells during creation of cohesive cells."
                << "  classifySize: " << classifySize << ", classifyTotal: " << classifyTotal;
            throw std::logic_error(msg.str());
        } // if
    }
    replaceCells.insert(vReplaceCells.begin(), vReplaceCells.end());
    // More checking
    noReplaceCells.insert(vNoReplaceCells.begin(), vNoReplaceCells.end());
} // classifyCellsDM
#endif

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

#if 0
    const PetscBool hasLagrangeConstraints = PETSC_TRUE;
    err = DMPlexCreateCohesiveSubmesh(dmDomain, hasLagrangeConstraints, labelName, labelValue, &dmFaultMesh));
#else

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
#endif
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
            PYLITH_JOURNAL_LOGICERROR("Inconsistent topology. Expected support of size 2 for face "
                                      << cone[iCone] << " of cohesive cell " <<cohesiveCell
                                      <<". Support has size "<<supportSize<<".");
        } // if
        assert(2 == supportSize);
        if ((cohesiveCell != support[0]) && (cohesiveCell != support[1]) ) {
            PYLITH_JOURNAL_LOGICERROR("Inconsistent topology. Cohesive cell "
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
