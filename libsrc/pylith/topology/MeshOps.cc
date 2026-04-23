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

#include "pylith/topology/MeshOps.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/utils/array.hh" // USES int_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/scales/Scales.hh" // USES Scales

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#include <algorithm> // USES std::sort, std::find
#include <map> // USES std::map

namespace pylith {
    namespace topology {
        class _MeshOps {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt createSubdomainMesh;
                static PylithInt createLowerDimMesh;
                static PylithInt createFromPoints;
                static PylithInt nondimensionalize;
                static PylithInt checkTopology;
                static PylithInt checkTopologyGeometry;
                static PylithInt checkTopologySymmetry;
                static PylithInt checkTopologySkeleton;
                static PylithInt checkMaterialLabels;

                static bool isInitialized;
            };
        };
    }
}
pylith::utils::EventLogger pylith::topology::_MeshOps::Events::logger;
PylithInt pylith::topology::_MeshOps::Events::createSubdomainMesh;
PylithInt pylith::topology::_MeshOps::Events::createLowerDimMesh;
PylithInt pylith::topology::_MeshOps::Events::createFromPoints;
PylithInt pylith::topology::_MeshOps::Events::nondimensionalize;
PylithInt pylith::topology::_MeshOps::Events::checkTopology;
PylithInt pylith::topology::_MeshOps::Events::checkTopologyGeometry;
PylithInt pylith::topology::_MeshOps::Events::checkTopologySymmetry;
PylithInt pylith::topology::_MeshOps::Events::checkTopologySkeleton;
PylithInt pylith::topology::_MeshOps::Events::checkMaterialLabels;
bool pylith::topology::_MeshOps::Events::isInitialized = false;

void
pylith::topology::_MeshOps::Events::init(void) {
    if (isInitialized) {
        return;
    } // if

    logger.setClassName("MeshOps");
    logger.initialize();
    createSubdomainMesh = logger.registerEvent("PL:MeshOps:createSubdomainMesh");
    createLowerDimMesh = logger.registerEvent("PL:MeshOps:createLowerDimMesh");
    createFromPoints = logger.registerEvent("PL:MeshOps:createFromPoints");
    nondimensionalize = logger.registerEvent("PL:MeshOps:nondimensionalize");
    checkTopology = logger.registerEvent("PL:MeshOps:checkTopology");
    checkTopologyGeometry = logger.registerEvent("PL:MeshOps:checkTopologyGeometry");
    checkTopologySymmetry = logger.registerEvent("PL:MeshOps:checkTopologySymmetry");
    checkTopologySkeleton = logger.registerEvent("PL:MeshOps:checkTopologySkeleton");
    checkMaterialLabels = logger.registerEvent("PL:MeshOps:checkMaterialLabels");

    isInitialized = true;
}


// ------------------------------------------------------------------------------------------------
// Create subdomain mesh using label.
pylith::topology::Mesh*
pylith::topology::MeshOps::createSubdomainMesh(const pylith::topology::Mesh& mesh,
                                               const char* labelName,
                                               const int labelValue,
                                               const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createSubdomainMesh);

    assert(labelName);

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);

    PetscBool hasLabel = PETSC_FALSE;
    PylithCallPetsc(DMHasLabel(dmDomain, labelName, &hasLabel));
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* :TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = NULL;
    PylithCallPetsc(DMGetLabel(dmDomain, labelName, &dmLabel));assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    PylithCallPetsc(DMLabelHasValue(dmLabel, labelValue, &hasLabelValue));
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    PylithCallPetsc(MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmDomain)));
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDM dmSubdomain = NULL;
    PylithCallPetsc(DMPlexFilter(dmDomain, dmLabel, labelValue, PETSC_FALSE, PETSC_FALSE, NULL, &dmSubdomain));

    PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
    PylithCallPetsc(DMPlexGetMaxSizes(dmSubdomain, &maxConeSizeLocal, NULL));
    PylithCallPetsc(MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmSubdomain)));

    if (maxConeSize <= 0) {
        PylithCallPetsc(DMDestroy(&dmSubdomain));
        std::ostringstream msg;
        msg << "Error while creating mesh of subdomain. Subdomain mesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    PylithScalar lengthScale;
    PylithCallPetsc(DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(DMPlexSetScale(dmSubdomain, PETSC_UNIT_LENGTH, lengthScale));

    pylith::topology::Mesh* submesh = new pylith::topology::Mesh();assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());

    submesh->setDM(dmSubdomain, componentName);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createSubdomainMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createSubdomainMesh


// ------------------------------------------------------------------------------------------------
// Create lower dimension mesh using label.
pylith::topology::Mesh*
pylith::topology::MeshOps::createLowerDimMesh(const pylith::topology::Mesh& mesh,
                                              const char* labelName,
                                              const int labelValue,
                                              const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createLowerDimMesh);
    assert(labelName);

    if (mesh.getDimension() < 1) {
        PYLITH_JOURNAL_LOGICERROR("Cannot create submesh for mesh with dimension < 1.");
    } // if

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscBool hasLabel = PETSC_FALSE;
    PylithCallPetsc(DMHasLabel(dmDomain, labelName, &hasLabel));
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = NULL;
    PylithCallPetsc(DMGetLabel(dmDomain, labelName, &dmLabel));assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    PylithCallPetsc(DMLabelHasValue(dmLabel, labelValue, &hasLabelValue));
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    PylithCallPetsc(MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmDomain)));
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    PetscInt labelHasVertices = 0;
    { // TEMPORARY: Continue to support creating lower dimension meshes using labels with vertices.
        PetscIS labelIS = NULL;
        const PetscInt* labelPoints = NULL;
        PetscInt numPoints = 0;
        PylithCallPetsc(DMGetStratumIS(dmDomain, labelName, labelValue, &labelIS));
        PetscInt labelHasVerticesLocal = 0;
        if (labelIS) {
            PylithCallPetsc(ISGetIndices(labelIS, &labelPoints));
            PylithCallPetsc(DMGetStratumSize(dmDomain, labelName, labelValue, &numPoints));

            topology::Stratum verticesStratum(dmDomain, topology::Stratum::DEPTH, 0);
            const PetscInt vStart = verticesStratum.begin();
            const PetscInt vEnd = verticesStratum.end();
            for (PetscInt iPoint = 0; iPoint < numPoints; ++iPoint) {
                if ((labelPoints[iPoint] >= vStart) && (labelPoints[iPoint] < vEnd) ) {
                    labelHasVerticesLocal = 1;
                    break;
                } // if
            } // if
            PylithCallPetsc(ISRestoreIndices(labelIS, &labelPoints));
        } // if
        PylithCallPetsc(ISDestroy(&labelIS));
        PylithCallPetsc(MPI_Allreduce(&labelHasVerticesLocal, &labelHasVertices, 1, MPI_INT, MPI_MAX,
                                      PetscObjectComm((PetscObject) dmDomain)));

        if (labelHasVertices) {
            pythia::journal::warning_t warning("deprecated");
            warning << pythia::journal::at(__HERE__)
                    << "DEPRECATION: Creating lower dimension mesh from label with vertices. "
                    << "This feature will be removed in v6.0. "
                    << "In the future, you will need to mark boundaries not vertices for boundary conditions."
                    << pythia::journal::endl;
        } // if
    } // TEMPORARY

    // We use DMPlexCreateSubmesh() instead of DMPlexFilter, because we want the submesh to have
    // domain cells hanging off of it, which allows us to project from the submesh to the domain mesh
    // to set boundary conditions using the auxiliary fields defined over the submesh.
    // DMPlexCreateSubmesh() requires a completed label.
    PetscDMLabel dmLabelFull = NULL;
    PylithCallPetsc(DMLabelDuplicate(dmLabel, &dmLabelFull));
    PylithCallPetsc(DMPlexLabelComplete(dmDomain, dmLabelFull));

    PetscDM dmSubmesh = NULL;
    const PetscBool markedFaces = !labelHasVertices ? PETSC_TRUE : PETSC_FALSE;
    PylithCallPetsc(DMPlexCreateSubmesh(dmDomain, dmLabelFull, labelValue, markedFaces, &dmSubmesh));
    PylithCallPetsc(DMLabelDestroy(&dmLabelFull));

    PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
    PylithCallPetsc(DMPlexGetMaxSizes(dmSubmesh, &maxConeSizeLocal, NULL));
    PylithCallPetsc(MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                                  PetscObjectComm((PetscObject) dmSubmesh)));

    if (maxConeSize <= 0) {
        PylithCallPetsc(DMDestroy(&dmSubmesh));
        std::ostringstream msg;
        msg << "Error while creating lower dimension mesh. Submesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Set length scale
    PylithScalar lengthScale;
    PylithCallPetsc(DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(DMPlexSetScale(dmSubmesh, PETSC_UNIT_LENGTH, lengthScale));
    pylith::topology::Mesh* submesh = new pylith::topology::Mesh(true);assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());

    submesh->setDM(dmSubmesh, componentName);

    // Check topology
    MeshOps::checkTopology(*submesh);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createLowerDimMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createLowerDimMesh


// ------------------------------------------------------------------------------------------------
// Create 0-dimension mesh from points.
pylith::topology::Mesh*
pylith::topology::MeshOps::createFromPoints(const PylithReal* points,
                                            const size_t numPoints,
                                            const spatialdata::geocoords::CoordSys* cs,
                                            const PylithReal lengthScale,
                                            MPI_Comm comm,
                                            const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createFromPoints);
    assert(cs);

    const int meshDim = 0;
    pylith::topology::Mesh* mesh = new pylith::topology::Mesh(meshDim, comm);assert(mesh);

    PetscDM dmPoints = NULL;
    const PetscInt depth = 0;
    PetscInt dmNumPoints[1];
    dmNumPoints[0] = numPoints;
    pylith::int_array dmConeSizes(0, numPoints);
    pylith::int_array dmCones(0, numPoints);
    pylith::int_array dmConeOrientations(0, numPoints);

    const size_t spaceDim = cs->getSpaceDim();

    PylithCallPetsc(DMPlexCreate(comm, &dmPoints));
    PylithCallPetsc(DMSetDimension(dmPoints, 0));
    PylithCallPetsc(DMSetCoordinateDim(dmPoints, spaceDim));
    if (numPoints > 0) {
        PylithCallPetsc(DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &dmConeSizes[0], &dmCones[0],
                                            &dmConeOrientations[0], points));
    } else {
        PetscInt empty[1];
        empty[0] = 0;
        PylithCallPetsc(DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &empty[0], &empty[0],
                                            &empty[0], points));
    } // if/else

    PetscSF sf = NULL;
    PylithCallPetsc(DMGetPointSF(dmPoints, &sf));
    PylithCallPetsc(PetscSFSetGraph(sf, numPoints, 0, NULL, PETSC_COPY_VALUES, NULL, PETSC_COPY_VALUES));

    mesh->setDM(dmPoints, componentName);

    mesh->setCoordSys(cs);

    PylithCallPetsc(DMPlexSetScale(mesh->getDM(), PETSC_UNIT_LENGTH, lengthScale));

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createFromPoints);
    PYLITH_METHOD_RETURN(mesh);
} // createFromPoints


// ------------------------------------------------------------------------------------------------
// Nondimensionalize the finite-element mesh.
void
pylith::topology::MeshOps::nondimensionalize(Mesh* const mesh,
                                             const pylith::scales::Scales& scales) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::nondimensionalize);

    assert(mesh);

    PetscVec coordVec = NULL;
    const PylithScalar lengthScale = scales.getLengthScale();

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    PylithCallPetsc(DMGetCoordinatesLocal(dmMesh, &coordVec));assert(coordVec);
    PylithCallPetsc(VecScale(coordVec, 1.0/lengthScale));
    PylithCallPetsc(DMPlexSetScale(dmMesh, PETSC_UNIT_LENGTH, lengthScale));
    PylithCallPetsc(DMViewFromOptions(dmMesh, NULL, "-pylith_nondim_dm_view"));

    const PetscInt dim = mesh->getDimension();
    if (dim >= 1) {
        const PylithReal avgDomainDim = computeAvgDomainDim(*mesh);
        const PylithReal avgDimTolerance = 0.01;
        if (avgDomainDim < avgDimTolerance) {
            std::ostringstream msg;
            msg << "Average domain dimension (" << avgDomainDim << ") is less than minimum tolerance ("
                << avgDimTolerance << "). This usually means the length scale (" << lengthScale << ") used in the "
                << "nondimensionalization needs to be smaller. The length scale should be on the order of the size "
                << "of the features controlling the displacement variations (fault length or domain size).";
            throw std::runtime_error(msg.str());
        } // if/else
    } // if

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::nondimensionalize);
    PYLITH_METHOD_END;
} // nondimensionalize


// ------------------------------------------------------------------------------------------------
// Strip out "ghost" cells hanging off mesh
PetscDM
pylith::topology::MeshOps::removeHangingCells(const PetscDM& dmMesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmClean = PETSC_NULLPTR;

    MPI_Comm comm = PetscObjectComm((PetscObject) dmMesh);
    pylith::topology::Stratum cells(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    PetscInt hasHangingCellsLocal = 0;
    if (cells.begin() != cells.end()) {
        DMPolytopeType cellType;
        PylithCallPetsc(DMPlexGetCellType(dmMesh, cells.begin(), &cellType));
        hasHangingCellsLocal = DMPolytopeTypeGetDim(cellType) < 0; // Hanging cells have dim == -1
    } // if
    PetscInt hasHangingCellsGlobal = 0;
    PylithCallPetsc(MPI_Allreduce(&hasHangingCellsLocal, &hasHangingCellsGlobal, 1, MPIU_INT, MPI_MAX, comm));
    if (hasHangingCellsGlobal) {
        // Create label over cells 1 dimension lower
        PetscDMLabel labelInclude = PETSC_NULLPTR;
        const PetscInt labelValue = 1;
        PylithCallPetsc(DMLabelCreate(comm, "no_hanging_cells", &labelInclude));
        pylith::topology::Stratum faces(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
        for (PetscInt face = faces.begin(); face < faces.end(); ++face) {
            PylithCallPetsc(DMLabelSetValue(labelInclude, face, labelValue));
        } // for

        PylithCallPetsc(DMPlexFilter(dmMesh, labelInclude, labelValue, PETSC_FALSE, PETSC_FALSE, PETSC_NULLPTR, &dmClean));
        PylithCallPetsc(DMLabelDestroy(&labelInclude));

        // Create section using subpoint map to ensure sections are consistent.
        PetscIS subpointIS = PETSC_NULLPTR;
        PetscSection sectionOld = PETSC_NULLPTR, sectionNew = PETSC_NULLPTR;
        PylithCallPetsc(DMPlexGetSubpointIS(dmClean, &subpointIS));
        PylithCallPetsc(DMGetLocalSection(dmMesh, &sectionOld));
        PylithCallPetsc(PetscSectionCreateSubmeshSection(sectionOld, subpointIS, &sectionNew));
        PylithCallPetsc(DMSetLocalSection(dmClean, sectionNew));
        PylithCallPetsc(PetscSectionDestroy(&sectionNew));

        PetscReal lengthScale = 0.0;
        PylithCallPetsc(DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale));
        PylithCallPetsc(DMPlexSetScale(dmClean, PETSC_UNIT_LENGTH, lengthScale));
    } else {
        dmClean = dmMesh;
        PylithCallPetsc(PetscObjectReference((PetscObject) dmClean));
    } // if/else

    PYLITH_METHOD_RETURN(dmClean);
}


// ------------------------------------------------------------------------------------------------
// Check topology of mesh.
void
pylith::topology::MeshOps::checkTopology(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopology);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    DMLabel subpointMap;
    PylithCallPetsc(DMPlexGetSubpointMap(dmMesh, &subpointMap));
    PetscInt cellHeight = subpointMap ? 1 : 0;

    PylithCallPetsc(DMViewFromOptions(dmMesh, NULL, "-pylith_checktopo_dm_view"));

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologyGeometry);
    PylithCallPetsc(DMPlexCheckGeometry(dmMesh));
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologyGeometry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySymmetry);
    PylithCallPetsc(DMPlexCheckSymmetry(dmMesh));
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologySymmetry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySkeleton);
    PylithCallPetsc(DMPlexCheckSkeleton(dmMesh, cellHeight));
    PylithCallPetsc(DMPlexCheckOrphanVertices(dmMesh));
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologySkeleton);

    /* Other check functions that we are not using:
     *
     * DMPlexCheckFaces() - not compatible with cohesive cells.
     *
     * DMPlexCheckInterfaceCones() - very slow
     */

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopology);
    PYLITH_METHOD_END;
} // checkTopology


// ------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isSimplexMesh(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    bool isSimplex = false;

    const PetscDM dm = mesh.getDM();
    PetscInt vStart = 0, vEnd = 0;
    PylithCallPetsc(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));
    if (vStart != vEnd) { // Test for simplex only works if we have points.
        PetscInt closureSize = 0;
        PetscInt* closure = NULL;
        const int dim = mesh.getDimension();

        PylithCallPetsc(DMPlexGetTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure));
        PetscInt numVertices = 0;
        for (PetscInt c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                ++numVertices;
            } // if
        } // for
        if (numVertices == dim+1) {
            isSimplex = PETSC_TRUE;
        } // if
        PylithCallPetsc(DMPlexRestoreTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure));
    } // if

    // Communicate result of isSimplex to all processes.
    int intSimplexLocal = isSimplex ? 1 : 0;
    int intSimplexGlobal = 0;
    MPI_Allreduce(&intSimplexLocal, &intSimplexGlobal, 1, MPI_INT, MPI_LOR, mesh.getComm());
    isSimplex = intSimplexGlobal == 1;

    PYLITH_METHOD_RETURN(isSimplex);
} // isSimplexMesh


// ------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isCohesiveCell(const PetscDM dm,
                                          const PetscInt cell) {
    bool isCohesive = false;

    DMPolytopeType ct;
    PylithCallPetsc(DMPlexGetCellType(dm, cell, &ct));
    if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_TRI_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_QUAD_PRISM_TENSOR)) { isCohesive = true; }

    return isCohesive;
} // isCohesiveCell


// ------------------------------------------------------------------------------------------------
// Get number of vertices in mesh.
PylithInt
pylith::topology::MeshOps::getNumVertices(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PylithInt nvertices = 0;
    PylithInt begin = 0, end = 0;
    PylithCallPetsc(DMPlexGetDepthStratum(dmMesh, 0, &begin, &end));
    nvertices = end-begin;

    PYLITH_METHOD_RETURN(nvertices);
}


// ------------------------------------------------------------------------------------------------
// Get number of cells in mesh.
PylithInt
pylith::topology::MeshOps::getNumCells(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscInt ncells = 0;
    PylithInt begin = 0, end = 0;
    const int cellHeight = 0;
    PylithCallPetsc(DMPlexGetHeightStratum(dmMesh, cellHeight, &begin, &end));
    ncells = end-begin;

    PYLITH_METHOD_RETURN(ncells);
}


// ------------------------------------------------------------------------------------------------
// Get number of vertices in a cell.
PylithInt
pylith::topology::MeshOps::getNumCorners(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscInt numCorners = 0;
    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    PetscInt cStart, cEnd, vStart, vEnd, closureSize, *closure = NULL;
    const int cellHeight = 0;
    PylithCallPetsc(DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd));
    PylithCallPetsc(DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd));
    if (cEnd > cStart) {
        PylithCallPetsc(DMPlexGetTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure));
        for (PetscInt c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {++numCorners;}
        } // for
        PylithCallPetsc(DMPlexRestoreTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure));
    } // if

    PYLITH_METHOD_RETURN(numCorners);
}


// ------------------------------------------------------------------------------------------------
/** Compute nominal dimension of domain based on mesh bounding box.
 *
 * @param[in] mesh Finite-element mesh.
 * @returns Average cell size.
 */
PylithReal
pylith::topology::MeshOps::computeAvgDomainDim(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    const PetscInt dim = mesh.getDimension();
    if (dim < 1) {
        PYLITH_METHOD_RETURN(0.0);
    } // if
    PylithReal coordMin[3];
    PylithReal coordMax[3];

    PylithCallPetsc(DMGetBoundingBox(mesh.getDM(), coordMin, coordMax));
    PylithReal volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= coordMax[i] - coordMin[i];
    } // for
    assert(dim > 0);
    const PylithReal avgDomainDim = pow(volume, 1.0/dim);

    PYLITH_METHOD_RETURN(avgDomainDim);
} // computeAvgDomainDim


// ------------------------------------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialLabels(const pylith::topology::Mesh& mesh,
                                               pylith::int_array& labelValues) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkMaterialLabels);

    // Create map with indices for each material
    const size_t numIds = labelValues.size();
    std::map<int, int> materialIndex;
    for (size_t i = 0; i < numIds; ++i) {
        materialIndex[labelValues[i]] = i;
    } // for

    int_array matCellCounts(numIds);
    matCellCounts = 0;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    Stratum cellsStratum(dmMesh, Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    PetscDMLabel materialsLabel = NULL;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    PylithCallPetsc(DMGetLabel(dmMesh, labelName, &materialsLabel));assert(materialsLabel);

    int *matBegin = &labelValues[0];
    int *matEnd = &labelValues[0] + labelValues.size();
    std::sort(matBegin, matEnd);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt matId;

        PylithCallPetsc(DMLabelGetValue(materialsLabel, c, &matId));
        if (matId < 0) {
            // :KLUDGE: Skip cells that are probably hybrid cells in halo
            // around fault that we currently ignore when looping over
            // materials (including cohesive cells).
            continue;
        } // if
        const int *result = std::find(matBegin, matEnd, matId);
        if (result == matEnd) {
            std::ostringstream msg;
            msg << "Material label_value '" << matId << "' for cell '" << c
                << "' does not match the label_value of any materials or interfaces.";
            throw std::runtime_error(msg.str());
        } // if

        const size_t matIndex = materialIndex[matId];
        assert(0 <= matIndex && matIndex < numIds);
        ++matCellCounts[matIndex];
    } // for

    // Make sure each material has cells.
    int_array matCellCountsAll(matCellCounts.size());
    PylithCallPetsc(MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
                                  matCellCounts.size(), MPI_INT, MPI_SUM, mesh.getComm()));
    for (size_t i = 0; i < numIds; ++i) {
        const int matId = labelValues[i];
        const size_t matIndex = materialIndex[matId];
        assert(0 <= matIndex && matIndex < numIds);
        if (matCellCountsAll[matIndex] <= 0) {
            std::ostringstream msg;
            msg << "No cells associated with material with id '" << matId << "'.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkMaterialLabels);
    PYLITH_METHOD_END;
} // checkMaterialIds


// End of file
