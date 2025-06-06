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
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <iostream> // USES std::cout
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::createFault(pylith::topology::Mesh* faultMesh,
                                         const pylith::topology::Mesh& mesh,
                                         PetscDMLabel surfaceLabel,
                                         const int surfaceLabelValue) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);
    PetscErrorCode err;

    faultMesh->setCoordSys(mesh.getCoordSys());
    PetscDM dmDomain = mesh.getDM();assert(dmDomain);

    // Convert fault to a DM
    PetscDM dmFault = PETSC_NULLPTR;
    const char *groupName = "";

    if (surfaceLabel) {
        err = PetscObjectGetName((PetscObject)surfaceLabel, &groupName);PYLITH_CHECK_ERROR(err);
    }

    PetscInt labelHasVertices = 0;
    { // TEMPORARY: Continue to support creating lower dimension meshes using labels with vertices.
        PetscIS labelIS = NULL;
        PetscInt labelHasVerticesLocal = 0;
        if (surfaceLabel) {
            const PetscInt* labelPoints = NULL;
            PetscInt numPoints = 0;
            err = DMGetStratumIS(dmDomain, groupName, surfaceLabelValue, &labelIS);PYLITH_CHECK_ERROR(err);assert(labelIS);
            err = ISGetIndices(labelIS, &labelPoints);PYLITH_CHECK_ERROR(err);
            err = DMGetStratumSize(dmDomain, groupName, surfaceLabelValue, &numPoints);PYLITH_CHECK_ERROR(err);

            pylith::topology::Stratum verticesStratum(dmDomain, pylith::topology::Stratum::DEPTH, 0);
            PetscInt vStart = 0, vEnd = 0;
            vStart = verticesStratum.begin();
            vEnd = verticesStratum.end();
            for (PetscInt iPoint = 0; iPoint < numPoints; ++iPoint) {
                if ((labelPoints[iPoint] >= vStart) && (labelPoints[iPoint] < vEnd) ) {
                    labelHasVerticesLocal = 1;
                    break;
                } // if
            } // if
            err = ISRestoreIndices(labelIS, &labelPoints);PYLITH_CHECK_ERROR(err);
        } // if
        err = ISDestroy(&labelIS);PYLITH_CHECK_ERROR(err);

        err = MPI_Allreduce(&labelHasVerticesLocal, &labelHasVertices, 1, MPI_INT, MPI_MAX,
                            PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);

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
        err = DMLabelDuplicate(surfaceLabel, &surfaceLabelFull);PYLITH_CHECK_ERROR(err);
        err = DMPlexLabelComplete(dmDomain, surfaceLabelFull);PYLITH_CHECK_ERROR(err);
    } // if

    const PetscBool markedFaces = !labelHasVertices ? PETSC_TRUE : PETSC_FALSE;
    err = DMPlexCreateSubmesh(dmDomain, surfaceLabelFull, surfaceLabelValue, markedFaces, &dmFault);PYLITH_CHECK_ERROR(err);
    err = DMLabelDestroy(&surfaceLabelFull);PYLITH_CHECK_ERROR(err);

    PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
    err = DMPlexGetMaxSizes(dmFault, &maxConeSizeLocal, NULL);PYLITH_CHECK_ERROR(err);
    err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmFault));PYLITH_CHECK_ERROR(err);

    if (maxConeSize <= 0) {
        err = DMDestroy(&dmFault);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Error while creating fault. Fault " << groupName << " with label value "
            << surfaceLabelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Check that no cells have all vertices on the fault
    if (surfaceLabel) {
        PetscIS subpointIS;
        const PetscInt *dmpoints;
        PetscInt defaultValue, cStart, cEnd, vStart, vEnd;

        err = DMLabelGetDefaultValue(surfaceLabel, &defaultValue);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSubpointIS(dmFault, &subpointIS);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetHeightStratum(dmFault, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepthStratum(dmDomain, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
        for (PetscInt c = cStart; c < cEnd; ++c) {
            PetscBool invalidCell = PETSC_TRUE;
            PetscInt *closure = NULL;
            PetscInt closureSize;

            err = DMPlexGetTransitiveClosure(dmDomain, dmpoints[c], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
                PetscInt value = 0;

                if ((closure[cl] < vStart) || (closure[cl] >= vEnd)) { continue;}
                err = DMLabelGetValue(surfaceLabel, closure[cl], &value);PYLITH_CHECK_ERROR(err);
                if (value == defaultValue) {invalidCell = PETSC_FALSE;break;}
            } // for
            err = DMPlexRestoreTransitiveClosure(dmDomain, dmpoints[c], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            if (invalidCell) {
                std::ostringstream msg;
                msg << "Ambiguous fault surface. Cell "<<dmpoints[c]<<" has all of its vertices on the fault.";
                err = ISRestoreIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
                err = DMDestroy(&dmFault);PYLITH_CHECK_ERROR(err);
                throw std::runtime_error(msg.str());
            } // if
        } // for
        err = ISRestoreIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
    } // if
    err = DMPlexOrient(dmFault);PYLITH_CHECK_ERROR(err);

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
    PetscInt dim, cMax, cStart, cEnd, numCohesiveCellsOld;
    PetscErrorCode err;

    // Have to remember the old number of cohesive cells
    err = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    cMax = cStart;
    for (PetscInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
        if (pylith::topology::MeshOps::isCohesiveCell(dm, cell)) { break; }
    } // for
    numCohesiveCellsOld = cEnd - cMax;
    // Create cohesive cells
    err = DMPlexGetSubpointMap(faultMesh.getDM(), &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelDuplicate(subpointMap, &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelClearStratum(label, mesh->getDimension());PYLITH_CHECK_ERROR(err);
    // Fix over-aggressive completion of boundary label
    err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);
    if (faultBdLabel && (dim > 2)) {
        PetscIS bdIS;
        const PetscInt *bd;
        PetscInt fStart, fEnd, n, i;

        err = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);PYLITH_CHECK_ERROR(err);
        err = DMLabelGetStratumIS(faultBdLabel, faultBdLabelValue, &bdIS);PYLITH_CHECK_ERROR(err);
        err = ISGetLocalSize(bdIS, &n);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(bdIS, &bd);PYLITH_CHECK_ERROR(err);
        for (i = 0; i < n; ++i) {
            const PetscInt p = bd[i];

            // Remove faces
            if ((p >= fStart) && (p < fEnd)) {
                const PetscInt *edges,   *verts, *supportA, *supportB;
                PetscInt numEdges, numVerts, supportSizeA, sA, supportSizeB, sB, val, bval, e, s;
                PetscBool found = PETSC_FALSE;

                err = DMLabelClearValue(faultBdLabel, p, faultBdLabelValue);PYLITH_CHECK_ERROR(err);
                // Remove the cross edge
                err = DMPlexGetCone(dm, p, &edges);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetConeSize(dm, p, &numEdges);PYLITH_CHECK_ERROR(err);
                if (numEdges != 3) {
                    std::ostringstream msg;
                    msg << "Internal error while creating fault mesh. Face "<<p<<" has "<<numEdges<<" edges != 3.";
                    throw std::logic_error(msg.str());
                }
                for (e = 0; e < numEdges; ++e) {
                    err = DMPlexGetCone(dm, edges[e], &verts);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetConeSize(dm, edges[e], &numVerts);PYLITH_CHECK_ERROR(err);
                    if (numVerts != 2) {
                        std::ostringstream msg;
                        msg << "Internal error while creating fault mesh. Edge "<<edges[e]<<" has "<<numVerts<<" vertices != 2.";
                        throw std::logic_error(msg.str());
                    }
                    err = DMPlexGetSupportSize(dm, verts[0], &supportSizeA);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetSupport(dm, verts[0], &supportA);PYLITH_CHECK_ERROR(err);
                    for (s = 0, sA = 0; s < supportSizeA; ++s) {
                        err = DMLabelGetValue(label, supportA[s], &val);PYLITH_CHECK_ERROR(err);
                        err = DMLabelGetValue(faultBdLabel, supportA[s], &bval);PYLITH_CHECK_ERROR(err);
                        if (( val >= 0) && ( bval >= 0) ) { ++sA;}
                    }
                    err = DMPlexGetSupportSize(dm, verts[1], &supportSizeB);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetSupport(dm, verts[1], &supportB);PYLITH_CHECK_ERROR(err);
                    for (s = 0, sB = 0; s < supportSizeB; ++s) {
                        err = DMLabelGetValue(label, supportB[s], &val);PYLITH_CHECK_ERROR(err);
                        err = DMLabelGetValue(faultBdLabel, supportB[s], &bval);PYLITH_CHECK_ERROR(err);
                        if (( val >= 0) && ( bval >= 0) ) { ++sB;}
                    }
                    if ((sA > 2) && (sB > 2)) {
                        err = DMLabelClearValue(faultBdLabel, edges[e], faultBdLabelValue);PYLITH_CHECK_ERROR(err);
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
        err = ISRestoreIndices(bdIS, &bd);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&bdIS);PYLITH_CHECK_ERROR(err);
    }
    // Completes the set of cells scheduled to be replaced
    err = DMPlexOrientLabel(dm, label);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelCohesiveComplete(dm, label, faultBdLabel, faultBdLabelValue, PETSC_FALSE, PETSC_FALSE, faultMesh.getDM());PYLITH_CHECK_ERROR(err);
    err = DMPlexConstructCohesiveCells(dm, label, NULL, &sdm);PYLITH_CHECK_ERROR(err);

    const char* interfaceLabelName = pylith::topology::Mesh::cells_label_name;
    err = DMGetLabel(sdm, interfaceLabelName, &mlabel);PYLITH_CHECK_ERROR(err);
    if (mlabel) {
        err = DMPlexGetHeightStratum(sdm, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        cMax = cStart;
        for (PetscInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
            if (pylith::topology::MeshOps::isCohesiveCell(sdm, cell)) { break; }
        }
        assert(cStart == cEnd || cEnd > cMax + numCohesiveCellsOld);
        for (PetscInt cell = cMax; cell < cEnd - numCohesiveCellsOld; ++cell) {
            PetscInt onBd;

            /* Eliminate hybrid cells on the boundary of the split from cohesive label,
             * they are marked with -(cell number) since the hybrid cell number aliases vertices in the old mesh */
            err = DMLabelGetValue(label, -cell, &onBd);PYLITH_CHECK_ERROR(err);
            // if (onBd == dim) continue;
            err = DMLabelSetValue(mlabel, cell, cohesiveLabelValue);PYLITH_CHECK_ERROR(err);
        }
    }
    err = DMLabelDestroy(&label);PYLITH_CHECK_ERROR(err);

    PetscReal lengthScale = 1.0;
    err = DMPlexGetScale(dm, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(sdm, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(sdm, NULL, "-pylith_cohesive_dm_view");PYLITH_CHECK_ERROR(err);
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
    PetscErrorCode err;

    faultMesh->setCoordSys(mesh.getCoordSys());

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscDM dmFaultMesh = NULL;

    const PetscBool hasLagrangeConstraints = PETSC_TRUE;
    err = DMPlexCreateCohesiveSubmesh(dmDomain, hasLagrangeConstraints, labelName, labelValue, &dmFaultMesh);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dmFaultMesh, NULL, "-pylith_fault_dm_view");PYLITH_CHECK_ERROR(err);
    err = DMPlexOrient(dmFaultMesh);PYLITH_CHECK_ERROR(err);
    std::string meshLabel = std::string("fault_") + std::string(surfaceLabel);

    PetscReal lengthScale = 1.0;
    err = DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmFaultMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    faultMesh->setDM(dmFaultMesh, meshLabel.c_str());

    PYLITH_METHOD_END;
} // createFaultParallel


// ------------------------------------------------------------------------------------------------
void
pylith::faults::TopologyOps::classifyCellsDM(PetscDM dmDomain,
                                             PetscInt vertex,
                                             const int depth,
                                             const int faceSize,
                                             PetscInt firstCohesiveCell,
                                             PointSet& replaceCells,
                                             PointSet& noReplaceCells,
                                             const int debug) {
    // Replace all cells on a given side of the fault with a vertex on the fault
    PointSet vReplaceCells;
    PointSet vNoReplaceCells;
    const PetscInt *support;
    PetscInt supportSize, s, classifyTotal = 0;
    PetscBool modified = PETSC_FALSE;
    PetscErrorCode err;

    if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
    err = DMPlexGetSupportSize(dmDomain, vertex, &supportSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSupport(dmDomain, vertex, &support);PYLITH_CHECK_ERROR(err);
    for (s = 0; s < supportSize; ++s) {
        const PetscInt point = support[s];

        if (point >= firstCohesiveCell) { return;}
        if (replaceCells.find(point)   != replaceCells.end()) {vReplaceCells.insert(point);}
        if (noReplaceCells.find(point) != noReplaceCells.end()) { vNoReplaceCells.insert(point);}
        modified = PETSC_TRUE;
        ++classifyTotal;
    }
    PetscInt classifySize = vReplaceCells.size() + vNoReplaceCells.size();

    while (modified && (classifySize < classifyTotal)) {
        modified = PETSC_FALSE;
        for (s = 0; s < supportSize; ++s) {
            const PetscInt point = support[s];
            PetscBool classified = PETSC_FALSE;

            if (debug) {
                const PetscInt *cone;
                PetscInt coneSize;

                std::cout << "Checking neighbor " << vertex << std::endl;
                err = DMPlexGetConeSize(dmDomain, vertex, &coneSize);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dmDomain, vertex, &cone);PYLITH_CHECK_ERROR(err);
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
                const PetscInt *coveringPoints;
                PetscInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                err = DMPlexGetMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                err = DMPlexRestoreMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
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
                const PetscInt *coveringPoints;
                PetscInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                err = DMPlexGetMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                err = DMPlexRestoreMeet(dmDomain, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
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
    PetscErrorCode err;
    PetscDMLabel interfacesLabel = NULL;

    const char* interfacesLabelName = TopologyOps::getInterfacesLabelName();
    PetscBool hasInterfacesLabel = PETSC_FALSE;
    if (DMHasLabel(dm, interfacesLabelName, &hasInterfacesLabel)) {
        err = DMGetLabel(dm, interfacesLabelName, &interfacesLabel);PYLITH_CHECK_ERROR(err);
    } else {
        PetscInt dim = 0;
        PetscInt pStart = 0;
        PetscInt pEnd = 0;
        PetscInt pMax = 0;

        err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);
        err = DMCreateLabel(dm, interfacesLabelName);PYLITH_CHECK_ERROR(err);
        err = DMGetLabel(dm, interfacesLabelName, &interfacesLabel);PYLITH_CHECK_ERROR(err);
        for (PylithInt iDim = 0; iDim <= dim; ++iDim) {
            err = DMPlexGetHeightStratum(dm, iDim, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
            err = DMPlexGetSimplexOrBoxCells(dm, iDim, NULL, &pMax);PYLITH_CHECK_ERROR(err);
            for (PylithInt p = pMax; p < pEnd; ++p) {
                err = DMLabelSetValue(interfacesLabel, p, 1);PYLITH_CHECK_ERROR(err);
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

    PetscErrorCode err = 0;

    const PylithInt* cone = NULL;
    err = DMPlexGetCone(dmDomain, cohesiveCell, &cone);PYLITH_CHECK_ERROR(err);
    PylithInt adjacentCells[2];
    for (PylithInt iCone = 0; iCone < 2; ++iCone) {
        const PylithInt* support = NULL;
        PylithInt supportSize = 0;

        err = DMPlexGetSupport(dmDomain, cone[iCone], &support);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSupportSize(dmDomain, cone[iCone], &supportSize);PYLITH_CHECK_ERROR(err);
        if (2 != supportSize) {
            PYLITH_JOURNAL_LOGICERROR("Inconsistent topology. Expected support of size 2 for cohesive cell "
                                      <<cohesiveCell<<". Support has size "<<supportSize<<".");
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
