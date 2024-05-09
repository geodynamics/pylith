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

#include "pylith/meshio/MeshBuilder.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscdmlabel.h" // USES PetscDMLabel

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace meshio {
        class _MeshBuilder {
public:

            static
            void faceFromCellSide(PetscInt* face,
                                  const PetscInt cell,
                                  const PetscInt side,
                                  PetscDM dmMesh);

            static
            void cellSideFromFace(PetscInt* cell,
                                  PetscInt* side,
                                  const PetscInt face,
                                  PetscDM dmMesh);

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt buildMesh;
                static PylithInt setGroup;
                static PylithInt setGroupAddPoints;

                static bool isInitialized;
            };
        };
    }
}
pylith::utils::EventLogger pylith::meshio::_MeshBuilder::Events::logger;
PylithInt pylith::meshio::_MeshBuilder::Events::buildMesh;
PylithInt pylith::meshio::_MeshBuilder::Events::setGroup;
PylithInt pylith::meshio::_MeshBuilder::Events::setGroupAddPoints;
bool pylith::meshio::_MeshBuilder::Events::isInitialized = false;

void
pylith::meshio::_MeshBuilder::Events::init(void) {
    if (isInitialized) {
        return;
    } // if

    logger.setClassName("MeshBuilder");
    logger.initialize();
    buildMesh = logger.registerEvent("PL:MeshBuilder:buildMesh");
    setGroup = logger.registerEvent("PL:MeshBuilder:setGroup");
    setGroupAddPoints = logger.registerEvent("PL:MeshBuilder:setGroupAddPoints");
    isInitialized = true;
}


// ------------------------------------------------------------------------------------------------
// Set vertices and cells in mesh.
void
pylith::meshio::MeshBuilder::buildMesh(topology::Mesh* mesh,
                                       scalar_array* coordinates,
                                       const int numVertices,
                                       int spaceDim,
                                       const int_array& cells,
                                       const int numCells,
                                       const int numCorners,
                                       const int meshDim,
                                       const bool isParallel) {
    PYLITH_METHOD_BEGIN;
    _MeshBuilder::Events::init();
    _MeshBuilder::Events::logger.eventBegin(_MeshBuilder::Events::buildMesh);

    assert(mesh);
    assert(coordinates);
    MPI_Comm comm = mesh->getComm();
    PetscInt dim = meshDim;
    PetscErrorCode err;

    { // Check to make sure every vertex is in at least one cell.
      // This is required by PETSc
        std::vector<bool> vertexInCell(numVertices, false);
        const int size = cells.size();
        for (int i = 0; i < size; ++i) {
            vertexInCell[cells[i]] = true;
        }
        int count = 0;
        for (int i = 0; i < numVertices; ++i) {
            if (!vertexInCell[i]) {
                ++count;
            }
        }
        if (count > 0) {
            std::ostringstream msg;
            msg << "Mesh contains " << count << " vertices that are not in any cells.";
            throw std::runtime_error(msg.str());
        } // if
    } // check

    /* DMPlex */
    PetscDM dmMesh = NULL;
    PetscBool interpolate = PETSC_TRUE;

    err = MPI_Bcast(&dim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
    err = MPI_Bcast(&spaceDim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
    const PetscInt bound = numCells*numCorners;
    for (PetscInt coff = 0; coff < bound; coff += numCorners) {
        DMPolytopeType ct;

        if (dim < 3) { continue;}
        switch (numCorners) {
        case 4: ct = DM_POLYTOPE_TETRAHEDRON;break;
        case 6: ct = DM_POLYTOPE_TRI_PRISM;break;
        case 8: ct = DM_POLYTOPE_HEXAHEDRON;break;
        default: continue;
        }
        err = DMPlexInvertCell(ct, (int *) &cells[coff]);PYLITH_CHECK_ERROR(err);
    } // for
    err = DMPlexCreateFromCellListPetsc(comm, dim, numCells, numVertices, numCorners, interpolate, &cells[0], spaceDim, &(*coordinates)[0], &dmMesh);PYLITH_CHECK_ERROR(err);
    mesh->setDM(dmMesh, "domain");

    _MeshBuilder::Events::logger.eventEnd(_MeshBuilder::Events::buildMesh);
    PYLITH_METHOD_END;
} // buildMesh


// ------------------------------------------------------------------------------------------------
// Build a point group as an int section.
void
pylith::meshio::MeshBuilder::setGroup(pylith::topology::Mesh* mesh,
                                      const char* name,
                                      const GroupPtType groupType,
                                      const int_array& points) {
    PYLITH_METHOD_BEGIN;
    _MeshBuilder::Events::init();
    _MeshBuilder::Events::logger.eventBegin(_MeshBuilder::Events::setGroup);
    assert(mesh);

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    const PetscInt numPoints = points.size();
    DMLabel dmLabel = PETSC_NULLPTR;
    const PetscInt labelValue = 1;
    PetscErrorCode err = PETSC_SUCCESS;

    err = DMCreateLabel(dmMesh, name);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, name, &dmLabel);PYLITH_CHECK_ERROR(err);
    if (CELL == groupType) {
        for (PetscInt p = 0; p < numPoints; ++p) {
            err = DMLabelSetValue(dmLabel, points[p], labelValue);PYLITH_CHECK_ERROR(err);
        } // for
    } else if (FACE == groupType) {
        assert(points.size() % 2 == 0);
        const size_t numFaces = points.size() / 2;
        for (size_t iFace = 0; iFace < numFaces; ++iFace) {
            const PylithInt cell = points[2*iFace+0];
            const PylithInt side = points[2*iFace+1];
            PetscInt face = -1;
            _MeshBuilder::faceFromCellSide(&face, cell, side, dmMesh);
            err = DMLabelSetValue(dmLabel, face, labelValue);PYLITH_CHECK_ERROR(err);
        } // for
    } else if (VERTEX == groupType) {
        PetscInt cStart, cEnd, vStart, vEnd, numCells;

        err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
        numCells = cEnd - cStart;
        for (PetscInt p = 0; p < numPoints; ++p) {
            err = DMLabelSetValue(dmLabel, numCells+points[p], labelValue);PYLITH_CHECK_ERROR(err);
        } // for
          // Also add any non-cells which have all vertices marked
        _MeshBuilder::Events::logger.eventBegin(_MeshBuilder::Events::setGroupAddPoints);
        for (PetscInt p = 0; p < numPoints; ++p) {
            const PetscInt vertex = numCells+points[p];
            PetscInt      *star = NULL, starSize, s;

            err = DMPlexGetTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
            for (s = 0; s < starSize*2; s += 2) {
                const PetscInt point = star[s];
                PetscInt      *closure = NULL, closureSize, c, value;
                PetscBool marked = PETSC_TRUE;

                if ((point >= cStart) && (point < cEnd)) { continue;}
                err = DMPlexGetTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                for (c = 0; c < closureSize*2; c += 2) {
                    if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                        err = DMLabelGetValue(dmLabel, closure[c], &value);PYLITH_CHECK_ERROR(err);
                        if (value != 1) {marked = PETSC_FALSE;break;}
                    }
                }
                err = DMPlexRestoreTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                if (marked) {err = DMLabelSetValue(dmLabel, point, labelValue);PYLITH_CHECK_ERROR(err);}
            }
            err = DMPlexRestoreTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
        }
        _MeshBuilder::Events::logger.eventEnd(_MeshBuilder::Events::setGroupAddPoints);
    } // if/else

    _MeshBuilder::Events::logger.eventEnd(_MeshBuilder::Events::setGroup);
    PYLITH_METHOD_END;
} // setGroup


// ------------------------------------------------------------------------------------------------
// Get a point group as an array of points.
void
pylith::meshio::MeshBuilder::getGroup(GroupPtType* groupType,
                                      int_array* points,
                                      const pylith::topology::Mesh& mesh,
                                      const char* name) {
    PYLITH_METHOD_BEGIN;

    assert(points);
    assert(groupType);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    const PetscInt numCells = cellsStratum.size();

    topology::Stratum facesStratum(dmMesh, topology::Stratum::HEIGHT, 1);
    const PetscInt fStart = facesStratum.begin();
    const PetscInt fEnd = facesStratum.end();

    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    PetscIS groupIS = NULL;
    const PetscInt* groupIndices = NULL;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMGetStratumIS(dmMesh, name, 1, &groupIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);

    PetscInt totalSize;
    err = DMGetStratumSize(dmMesh, name, 1, &totalSize);PYLITH_CHECK_ERROR(err);

    *groupType = pylith::meshio::MeshBuilder::VERTEX;
    if (( totalSize > 0) && (( groupIndices[0] >= cStart) && (groupIndices[0] < cEnd) )) {
        *groupType = pylith::meshio::MeshBuilder::CELL;
    } else if (( totalSize > 0) && (( groupIndices[0] >= fStart) && (groupIndices[0] < fEnd) )) {
        *groupType = pylith::meshio::MeshBuilder::FACE;
    } // if

    PetscInt offset = 0;
    PetscInt pStart = 0;
    PetscInt pEnd = 0;
    switch (*groupType) {
    case pylith::meshio::MeshBuilder::VERTEX:
        offset = numCells;
        pStart = vStart;
        pEnd = vEnd;
        break;
    case pylith::meshio::MeshBuilder::FACE:
        offset = 0;
        pStart = fStart;
        pEnd = fEnd;
        break;
    case pylith::meshio::MeshBuilder::CELL:
        offset = 0;
        pStart = cStart;
        pEnd = cEnd;
        break;
    } // switch

    // Count number of vertices/faces/cells, filtering out any points not in group type.
    PetscInt groupSize = 0;
    for (PetscInt p = 0; p < totalSize; ++p) {
        if (( groupIndices[p] >= pStart) && ( groupIndices[p] < pEnd) ) {
            ++groupSize;
        } // if
    } // for

    points->resize(groupSize);
    for (PetscInt p = 0; p < groupSize; ++p) {
        if (( groupIndices[p] >= pStart) && ( groupIndices[p] < pEnd) ) {
            (*points)[p] = groupIndices[p]-offset;
        } // if
    } // for
    err = ISRestoreIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&groupIS);PYLITH_CHECK_ERROR(err);

    // Convert face point to cell+side.
    if (pylith::meshio::MeshBuilder::FACE == *groupType) {
        int_array faces(groupSize);
        faces = *points;
        points->resize(2*groupSize);
        for (int iFace = 0; iFace < groupSize; ++iFace) {
            _MeshBuilder::cellSideFromFace(&(*points)[2*iFace+0], &(*points)[2*iFace+1], faces[iFace], dmMesh);
        } // for
    } // if

    PYLITH_METHOD_END;
} // getGroup


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_MeshBuilder::faceFromCellSide(PetscInt* face,
                                               const PetscInt cell,
                                               const PetscInt side,
                                               PetscDM dmMesh) {
    PYLITH_METHOD_BEGIN;
    assert(face);

    PetscErrorCode err = PETSC_SUCCESS;

    DMPolytopeType cellType;
    err = DMPlexGetCellType(dmMesh, cell, &cellType);PYLITH_CHECK_ERROR(err);

    const PetscInt* cone = NULL;
    PetscInt coneSize = 0;
    err = DMPlexGetCone(dmMesh, cell, &cone);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(dmMesh, cell, &coneSize);PYLITH_CHECK_ERROR(err);
    PetscInt coneIndex = -1;
    switch (cellType) {
    case DM_POLYTOPE_TRIANGLE: {
        // side: 0=bottom, 1=right, 2=left
        // face: 0=bottom, 1=right, 2=left
        const int faceMapping[3] = { 0, 1, 2 };
        assert(coneSize == 3);
        coneIndex = faceMapping[side];
        break;
    } // triangle
    case DM_POLYTOPE_QUADRILATERAL: {
        // side: 0=bottom, 1=right, 2=top, 3=left
        // face: 0=bottom, 1=right, 2=top, 3=left
        const int faceMapping[4] = { 0, 1, 2, 3 };
        assert(coneSize == 4);
        coneIndex = faceMapping[side];
        break;
    } // quadrilateral
    case DM_POLYTOPE_TETRAHEDRON: {
        // side:  0  1  2  3
        // face:  1  3  2  0
        const int faceMapping[4] = { 1, 3, 2, 0 };
        assert(coneSize == 4);
        coneIndex = faceMapping[side];
        break;
    } // tetrahedron
    case DM_POLYTOPE_HEXAHEDRON: {
        // side: 0=front, 1=right, 2=back, 3=left, 4=bottom, 5=top
        // face: 2=front, 4=right, 3=back, 5=left, 0=bottom, 1=top
        const int faceMapping[6] = { 2, 4, 3, 5, 0, 1 };
        assert(coneSize == 6);
        coneIndex = faceMapping[side];
        break;
    } // hexahedron

    case DM_POLYTOPE_SEG_PRISM_TENSOR:
    case DM_POLYTOPE_POINT_PRISM_TENSOR:
    case DM_POLYTOPE_SEGMENT:
    case DM_POLYTOPE_POINT:
    default:
        PYLITH_JOURNAL_LOGICERROR("Can only get faces for triangles, quadrilaterals, tetrahedra, and hexahedra.");
    } // switch
    assert(coneIndex >= 0);
    *face = cone[coneIndex];

    PYLITH_METHOD_END;
} // faceFromCellSide


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_MeshBuilder::cellSideFromFace(PetscInt* cell,
                                               PetscInt* side,
                                               const PetscInt face,
                                               PetscDM dmMesh) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscInt* support = NULL;
    PetscInt supportSize = 0;
    err = DMPlexGetSupport(dmMesh, face, &support);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSupportSize(dmMesh, face, &supportSize);PYLITH_CHECK_ERROR(err);
    if (!supportSize || !support) {
        PYLITH_JOURNAL_LOGICERROR("Could not find support for face " << face << ".");
    } // if
    *cell = support[0];
    DMPolytopeType cellType;
    err = DMPlexGetCellType(dmMesh, *cell, &cellType);PYLITH_CHECK_ERROR(err);

    // Find side of cell using cone.
    const PetscInt* cone = NULL;
    PetscInt coneSize = 0;
    PetscInt coneIndex = -1;
    err = DMPlexGetCone(dmMesh, *cell, &cone);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetConeSize(dmMesh, *cell, &coneSize);PYLITH_CHECK_ERROR(err);
    for (PetscInt iCone = 0; iCone < coneSize; ++iCone) {
        if (cone[iCone] == face) {
            coneIndex = iCone;
            break;
        } // if
    } // for
    if (coneIndex < 0) {
        PYLITH_JOURNAL_LOGICERROR("Could not find face '" << face << "' in cone of cell '" << *cell << "'.");
    } // if

    switch (cellType) {
    case DM_POLYTOPE_TRIANGLE: {
        // face: 0=bottom, 1=right, 2=left
        // side: 0=bottom, 1=right, 2=left
        const int sideMapping[3] = { 0, 1, 2 };
        *side = sideMapping[coneIndex];
        break;
    } // triangle
    case DM_POLYTOPE_QUADRILATERAL: {
        // face: 0=bottom, 1=right, 2=top, 3=left
        // side: 0=bottom, 1=right, 2=top, 3=left
        const int sideMapping[4] = { 0, 1, 2, 3 };
        *side = sideMapping[coneIndex];
        break;
    } // quadrilateral
    case DM_POLYTOPE_TETRAHEDRON: {
        // face:  0  1  2  3
        // side:  3  0  2  1
        const int sideMapping[4] = { 3, 0, 2, 1 };
        *side = sideMapping[coneIndex];
        break;
    } // tetrahedron
    case DM_POLYTOPE_HEXAHEDRON: {
        // face: 0=bottom, 1=top, 2=front, 3=back, 4=right, 5=left
        // side: 4=bottom, 5=top, 0=front, 2=back, 1=right, 3=left
        const int sideMapping[6] = { 4, 5, 0, 2, 1, 3 };
        *side = sideMapping[coneIndex];
        break;
    } // hexahedron

    case DM_POLYTOPE_SEG_PRISM_TENSOR:
    case DM_POLYTOPE_POINT_PRISM_TENSOR:
    case DM_POLYTOPE_SEGMENT:
    case DM_POLYTOPE_POINT:
    default:
        PYLITH_JOURNAL_LOGICERROR("Can only get faces for triangles, quadrilaterals, tetrahedra, and hexahedra.");
    } // switch

    PYLITH_METHOD_END;
}


// End of file
