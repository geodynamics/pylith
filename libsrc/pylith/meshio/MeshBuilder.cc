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
#include "pylith/topology/MeshOps.hh" // USES MeshOps
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
            void getGroupNames(string_vector* names,
                               const pylith::topology::Mesh& mesh,
                               const PetscInt pStart,
                               const PetscInt pEnd,
                               const bool exclusive);

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

            static
            void faceFromCellVertices(PetscInt* face,
                                      const PetscInt cell,
                                      const PetscInt* faceVertices,
                                      const PetscInt numFaceVertices,
                                      PetscDM dmMesh);

            static
            void cellVerticesFromFace(PetscInt* cell,
                                      int_array* faceVertices,
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
                                       const Topology& topology,
                                       const Geometry& geometry,
                                       const bool isParallel) {
    PYLITH_METHOD_BEGIN;
    _MeshBuilder::Events::init();
    _MeshBuilder::Events::logger.eventBegin(_MeshBuilder::Events::buildMesh);

    assert(mesh);
    MPI_Comm comm = mesh->getComm();
    PetscInt dim = topology.dimension;
    PetscErrorCode err;

    { // Check to make sure every vertex is in at least one cell.
      // This is required by PETSc
        std::vector<bool> vertexInCell(geometry.numVertices, false);
        const size_t size = topology.cells.size();
        for (size_t i = 0; i < size; ++i) {
            vertexInCell[topology.cells[i]] = true;
        }
        size_t count = 0;
        for (size_t i = 0; i < geometry.numVertices; ++i) {
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

    err = MPI_Bcast(&dim, 1, MPIU_INT, 0, comm);PYLITH_CHECK_ERROR(err);
    const PetscInt bound = topology.numCells * topology.numCorners;
    int_array cellsCopy(topology.cells); // Use copy because we reuse in testing.
    if (3 == topology.dimension) {
        DMPolytopeType ct;
        switch (topology.cellShape) {
        case TETRAHEDRON: ct = DM_POLYTOPE_TETRAHEDRON;break;
        case HEXAHEDRON: ct = DM_POLYTOPE_HEXAHEDRON;break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown cell shape.");
        }
        for (PetscInt coff = 0; coff < bound; coff += topology.numCorners) {
            err = DMPlexInvertCell(ct, (int *) &cellsCopy[coff]);PYLITH_CHECK_ERROR(err);
        } // for
    } // if

    PetscDM dmMesh = NULL;
    PetscBool interpolate = PETSC_TRUE;
    err = DMPlexCreateFromCellListPetsc(comm, dim, topology.numCells, geometry.numVertices, topology.numCorners, interpolate, &cellsCopy[0], dim, &geometry.vertices[0], &dmMesh);PYLITH_CHECK_ERROR(err);
    mesh->setDM(dmMesh, "domain");

    _MeshBuilder::Events::logger.eventEnd(_MeshBuilder::Events::buildMesh);
    PYLITH_METHOD_END;
} // buildMesh


// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshBuilder::setMaterials(pylith::topology::Mesh* mesh,
                                          const int_array& materialIds) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    PetscErrorCode err = PETSC_SUCCESS;
    if (!mesh->getCommRank()) {
        pylith::topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();

        if (size_t(cellsStratum.size()) != materialIds.size()) {
            std::ostringstream msg;
            msg << "Mismatch in size of materials identifier array ("
                << materialIds.size() << ") and number of cells in mesh ("<< (cEnd - cStart) << ").";
            throw std::runtime_error(msg.str());
        } // if
        for (PetscInt c = cStart; c < cEnd; ++c) {
            err = DMSetLabelValue(dmMesh, labelName, c, materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
        } // for
    } else {
        err = DMCreateLabel(dmMesh, labelName);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setMaterials


// ------------------------------------------------------------------------------------------------
// Build a point group for vertices.
void
pylith::meshio::MeshBuilder::setVertexGroup(pylith::topology::Mesh* mesh,
                                            const char* name,
                                            const int_array& points,
                                            const int labelValue) {
    PYLITH_METHOD_BEGIN;
    _MeshBuilder::Events::init();
    _MeshBuilder::Events::logger.eventBegin(_MeshBuilder::Events::setGroup);
    assert(mesh);

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    const PetscInt numPoints = points.size();
    DMLabel dmLabel = PETSC_NULLPTR;
    PetscErrorCode err = PETSC_SUCCESS;

    err = DMCreateLabel(dmMesh, name);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, name, &dmLabel);PYLITH_CHECK_ERROR(err);

    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    const PetscInt offset = cellsStratum.size();

    pylith::topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    for (PetscInt p = 0; p < numPoints; ++p) {
        err = DMLabelSetValue(dmLabel, offset+points[p], labelValue);PYLITH_CHECK_ERROR(err);
    } // for
    // Also add any non-cells which have all vertices marked
    for (PetscInt p = 0; p < numPoints; ++p) {
        const PetscInt vertex = offset+points[p];
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
                } // if
            } // for
            err = DMPlexRestoreTransitiveClosure(dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            if (marked) {err = DMLabelSetValue(dmLabel, point, labelValue);PYLITH_CHECK_ERROR(err);}
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
    } // for

    _MeshBuilder::Events::logger.eventEnd(_MeshBuilder::Events::setGroup);
    PYLITH_METHOD_END;
} // setVertexGroup


// ------------------------------------------------------------------------------------------------
// Build a point group for faces from cell+vertices.
void
pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(pylith::topology::Mesh* mesh,
                                                          const char* name,
                                                          const int_array& faceValues,
                                                          const shape_t faceShape,
                                                          const int labelValue) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    const size_t numFaceVertices = getNumVerticesFace(faceShape);
    size_t numFaceValues = 1 + numFaceVertices;
    assert(0 == faceValues.size() % numFaceValues);
    const size_t numFaces = faceValues.size() / numFaceValues;

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    DMLabel dmLabel = PETSC_NULLPTR;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMCreateLabel(dmMesh, name);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, name, &dmLabel);PYLITH_CHECK_ERROR(err);

    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt offset = cellsStratum.size();

    for (size_t iFace = 0; iFace < numFaces; ++iFace) {
        const PylithInt cell = faceValues[numFaceValues*iFace+0];
        int_array faceVertices(&faceValues[numFaceValues*iFace+1], numFaceValues-1);
        faceVertices += offset;
        PetscInt face = -1;
        _MeshBuilder::faceFromCellVertices(&face, cell, &faceVertices[0], numFaceValues-1, dmMesh);
        err = DMLabelSetValue(dmLabel, face, labelValue);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // setFaceGroupFromCellVertices


// ------------------------------------------------------------------------------------------------
// Build a point group for faces from cell+side.
void
pylith::meshio::MeshBuilder::setFaceGroupFromCellSide(pylith::topology::Mesh* mesh,
                                                      const char* name,
                                                      const int_array& faceValues,
                                                      const int labelValue) {
    PYLITH_METHOD_BEGIN;
    assert(mesh);

    size_t numFaceValues = 2; // cell + side
    assert(0 == faceValues.size() % numFaceValues);
    const size_t numFaces = faceValues.size() / numFaceValues;

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    DMLabel dmLabel = PETSC_NULLPTR;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMCreateLabel(dmMesh, name);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, name, &dmLabel);PYLITH_CHECK_ERROR(err);

    for (size_t index = 0, iFace = 0; iFace < numFaces; ++iFace) {
        const PylithInt cell = faceValues[index++];
        const PylithInt side = faceValues[index++];
        PetscInt face = -1;
        _MeshBuilder::faceFromCellSide(&face, cell, side, dmMesh);
        err = DMLabelSetValue(dmLabel, face, labelValue);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // setFaceGroupFromCellSide


// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshBuilder::getVertices(Geometry* geometry,
                                         const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    assert(geometry);

    geometry->spaceDim = mesh.getDimension();

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscVec coordVec = NULL;
    PetscScalar* coordArray = NULL;
    PetscInt coordSize = 0;
    PylithScalar lengthScale = 1.0;
    PetscErrorCode err = 0;

    // Get length scale for dimensioning
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);

    // Get coordinates and dimensionalize values
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);
    err = VecGetLocalSize(coordVec, &coordSize);PYLITH_CHECK_ERROR(err);
    assert(0 == coordSize % geometry->spaceDim);
    geometry->numVertices = coordSize / geometry->spaceDim;

    geometry->vertices.resize(coordSize);
    for (PetscInt i = 0; i < coordSize; ++i) {
        geometry->vertices[i] = coordArray[i]*lengthScale;
    } // for
    err = VecRestoreArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _getVertices


// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshBuilder::getCells(Topology* topology,
                                      const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    assert(topology);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    topology->dimension = mesh.getDimension();
    topology->numCells = pylith::topology::MeshOps::getNumCells(mesh);assert(topology->numCells > 0);
    topology->numCorners = pylith::topology::MeshOps::getNumCorners(mesh);assert(topology->numCorners > 0);
    assert(size_t(cellsStratum.size()) == topology->numCells);
    topology->cellShape = cellShapeFromCorners(topology->dimension, topology->numCorners);

    topology->cells.resize(topology->numCells * topology->numCorners);

    PetscIS globalVertexNumbers = NULL;
    const PetscInt* gvertex = NULL;
    PetscErrorCode err = 0;

    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        DMPolytopeType ct;
        PetscInt numCorners = 0, closureSize, *closure = NULL;

        err = DMPlexGetCellType(dmMesh, c, &ct);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
            if ((closure[cl] >= vStart) && (closure[cl] < vEnd)) {
                const PetscInt gv = gvertex[closure[cl]-vStart];
                topology->cells[index++] = gv < 0 ? -(gv+1) : gv;
                ++numCorners;
            }
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        err = DMPlexInvertCell(ct, &topology->cells[index-numCorners]);PYLITH_CHECK_ERROR(err);
        assert(size_t(numCorners) == topology->numCorners);
    } // for
    err = ISRestoreIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // getCells


// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshBuilder::getMaterials(int_array* materialIds,
                                          const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    assert(materialIds);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    materialIds->resize(cellsStratum.size());
    PetscErrorCode err = PETSC_SUCCESS;
    PetscInt matId = 0;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, labelName, c, &matId);PYLITH_CHECK_ERROR(err);
        (*materialIds)[index++] = matId;
    } // for

    PYLITH_METHOD_END;
} // getMaterials


// ----------------------------------------------------------------------
// Get names of vertex groups in mesh.
void
pylith::meshio::MeshBuilder::getVertexGroupNames(string_vector* names,
                                                 const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    assert(names);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();
    const bool exclusive = false;
    _MeshBuilder::getGroupNames(names, mesh, vStart, vEnd, exclusive);

    PYLITH_METHOD_END;
} // getVertexGroupNames


// ----------------------------------------------------------------------
// Get names of face groups in mesh.
void
pylith::meshio::MeshBuilder::getFaceGroupNames(string_vector* names,
                                               const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    assert(names);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum facesStratum(dmMesh, topology::Stratum::HEIGHT, 1);
    const PetscInt fStart = facesStratum.begin();
    const PetscInt fEnd = facesStratum.end();
    const bool exclusive = true; // Ignore groups with points that are not faces
    _MeshBuilder::getGroupNames(names, mesh, fStart, fEnd, exclusive);

    PYLITH_METHOD_END;
} // getFaceGroupNames


// ------------------------------------------------------------------------------------------------
// Get a face group as an array of points.
void
pylith::meshio::MeshBuilder::getVertexGroup(int_array* points,
                                            const pylith::topology::Mesh& mesh,
                                            const char* name,
                                            const int labelValue) {
    PYLITH_METHOD_BEGIN;
    assert(points);
    assert(name);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt numCells = cellsStratum.size();

    pylith::topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    pylith::topology::StratumIS groupIS(dmMesh, name, labelValue);
    const PetscInt* groupPoints = groupIS.points();
    const PetscInt totalSize = groupIS.size();
    const PetscInt offset = numCells;

    PetscInt numPoints = 0;
    int_array buffer(totalSize);
    for (PetscInt p = 0; p < totalSize; ++p) {
        if (( groupPoints[p] >= vStart) && ( groupPoints[p] < vEnd) ) {
            buffer[numPoints++] = groupPoints[p]-offset;
        } // if
    } // for
    *points = int_array(&buffer[0], numPoints);

    PYLITH_METHOD_END;
} // getVertexGroup


// ------------------------------------------------------------------------------------------------
// Get a face group as an array of points.
void
pylith::meshio::MeshBuilder::getFaceGroup(int_array* points,
                                          const pylith::topology::Mesh& mesh,
                                          const char* name,
                                          const int labelValue) {
    PYLITH_METHOD_BEGIN;
    assert(points);
    assert(name);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt offset = cellsStratum.size();

    pylith::topology::Stratum facesStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
    const PetscInt fStart = facesStratum.begin();
    const PetscInt fEnd = facesStratum.end();

    pylith::topology::StratumIS groupIS(dmMesh, name, labelValue);
    const PetscInt* groupPoints = groupIS.points();
    const PetscInt totalSize = groupIS.size();

    PetscInt numFaces = 0;
    int_array buffer(totalSize);
    for (PetscInt p = 0; p < totalSize; ++p) {
        if (( groupPoints[p] >= fStart) && ( groupPoints[p] < fEnd) ) {
            buffer[numFaces++] = groupPoints[p];
        } // if
    } // for
    int_array faces(&buffer[0], numFaces);
    groupIS.deallocate();

    // Convert face point to cell+vertices.
    const PetscInt maxNumFaceVertices = 4;
    size_t index = 0;
    buffer.resize(numFaces*(1+maxNumFaceVertices));
    for (PetscInt iFace = 0; iFace < numFaces; ++iFace) {
        PetscInt cell;
        int_array faceVertices;
        _MeshBuilder::cellVerticesFromFace(&cell, &faceVertices, faces[iFace], dmMesh);
        buffer[index++] = cell;
        for (size_t i = 0; i < faceVertices.size(); ++i) {
            buffer[index++] = faceVertices[i] - offset;
        } // for
    } // for
    *points = int_array(&buffer[0], index);

    PYLITH_METHOD_END;
} // getFaceGroup


// ------------------------------------------------------------------------------------------------
// Get cell shape from dimension and number of corners.
pylith::meshio::MeshBuilder::shape_t
pylith::meshio::MeshBuilder::cellShapeFromCorners(const size_t cellDim,
                                                  const size_t numCorners) {
    shape_t cellShape = POINT;
    if ((2 == cellDim) && (3 == numCorners)) {
        cellShape = TRIANGLE;
    } else if ((2 == cellDim) && (4 == numCorners)) {
        cellShape = QUADRILATERAL;
    } else if ((3 == cellDim) && (4 == numCorners)) {
        cellShape = TETRAHEDRON;
    } else if ((3 == cellDim) && (8 == numCorners)) {
        cellShape = HEXAHEDRON;
    } else {
        PYLITH_JOURNAL_LOGICERROR("Unknown cell type. dim=" << cellDim <<", # corners="<<numCorners<<".");
    } // if/else

    return cellShape;
} // cellShapeFromCorners


// ------------------------------------------------------------------------------------------------
// Get face shape from cell shape.
pylith::meshio::MeshBuilder::shape_t
pylith::meshio::MeshBuilder::faceShapeFromCellShape(const shape_t cellShape) {
    shape_t faceShape = POINT;
    switch (cellShape) {
    case LINE: faceShape = POINT;break;
    case TRIANGLE: faceShape = LINE;break;
    case QUADRILATERAL: faceShape = LINE;break;
    case TETRAHEDRON: faceShape = TRIANGLE;break;
    case HEXAHEDRON: faceShape = QUADRILATERAL;break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown cell shape '"<<cellShape<<"'.");
    } // switch

    return faceShape;
} // faceShapeFromCellShape


// ------------------------------------------------------------------------------------------------
// Get number of face vertices given face shape.
size_t
pylith::meshio::MeshBuilder::getNumVerticesFace(const shape_t faceShape) {
    size_t numVertices = 0;
    switch (faceShape) {
    case POINT: numVertices = 1;break;
    case LINE: numVertices = 2;break;
    case TRIANGLE: numVertices = 3;break;
    case QUADRILATERAL: numVertices = 4;break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown face shape '" << faceShape << "'.");
    } // switch

    return numVertices;
}


// ----------------------------------------------------------------------
// Get names of groups in mesh with points in given range.
void
pylith::meshio::_MeshBuilder::getGroupNames(string_vector* names,
                                            const pylith::topology::Mesh& mesh,
                                            const PetscInt pStart,
                                            const PetscInt pEnd,
                                            const bool exclusive) {
    PYLITH_METHOD_BEGIN;
    assert(names);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscInt numLabels = 0;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMGetNumLabels(dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);

    const std::string& materialLabelName = pylith::topology::Mesh::cells_label_name;

    size_t numNames = 0;
    names->resize(numLabels);
    for (int iLabel = 0; iLabel < numLabels; ++iLabel) {
        const char* labelStr = NULL;
        err = DMGetLabelName(dmMesh, iLabel, &labelStr);PYLITH_CHECK_ERROR(err);
        const std::string labelName = std::string(labelStr);

        if ((labelName != std::string("depth"))
            && (labelName != std::string("celltype"))
            && (labelName != materialLabelName)) {
            PetscDMLabel dmLabel = PETSC_NULLPTR;
            err = DMGetLabel(dmMesh, labelStr, &dmLabel);PYLITH_CHECK_ERROR(err);
            PetscInt numLabelValues;
            PetscIS labelValuesIS = PETSC_NULLPTR;
            const PetscInt* labelValues = PETSC_NULLPTR;
            err = DMLabelGetNumValues(dmLabel, &numLabelValues);PYLITH_CHECK_ERROR(err); // assert(1 == numLabelValues);
            err = DMLabelGetValueIS(dmLabel, &labelValuesIS);PYLITH_CHECK_ERROR(err);assert(labelValuesIS);
            err = ISGetIndices(labelValuesIS, &labelValues);PYLITH_CHECK_ERROR(err);
            if (labelValues) {
                const PetscInt labelValue = labelValues[0];
                pylith::topology::StratumIS pointIS(dmMesh, labelStr, labelValue);
                const PetscInt* points = pointIS.points();
                const PetscInt numPoints = pointIS.size();
                bool hasOtherPoints = false;
                bool foundPoint = false;
                for (PetscInt iPoint = 0; iPoint < numPoints; ++iPoint) {
                    if ((points[iPoint] >= pStart) && (points[iPoint] < pEnd) ) {
                        foundPoint = true;
                        if (!exclusive) {
                            break;
                        } // if
                    } else {
                        hasOtherPoints = true;
                        if (exclusive) {
                            break;
                        } // if
                    } // if/else
                } // for
                if ((foundPoint && !exclusive) || (foundPoint && exclusive && !hasOtherPoints)) {
                    (*names)[numNames++] = labelName;
                } // if
            }
            err = ISRestoreIndices(labelValuesIS, &labelValues);PYLITH_CHECK_ERROR(err);
            err = ISDestroy(&labelValuesIS);PYLITH_CHECK_ERROR(err);
        } // if
    } // for
    names->resize(numNames);

    PYLITH_METHOD_END;
} // getGroupNames


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
    if ((side < 0) || (side >= coneSize)) {
        std::ostringstream msg;
        msg << "Cell side '" << side << "' must be in [0, " << coneSize << ").";
        throw std::runtime_error(msg.str());
    } // if
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
        // face:  1  2  3  0
        const int faceMapping[4] = { 1, 2, 3, 0 };
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
pylith::meshio::_MeshBuilder::faceFromCellVertices(PetscInt* face,
                                                   const PetscInt cell,
                                                   const PetscInt* faceVertices,
                                                   const PetscInt numFaceVertices,
                                                   PetscDM dmMesh) {
    PYLITH_METHOD_BEGIN;
    assert(face);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscInt numFaces = 0;
    const PetscInt* faces = PETSC_NULLPTR;
    err = DMPlexGetFullJoin(dmMesh, numFaceVertices, faceVertices, &numFaces, &faces);PYLITH_CHECK_ERROR(err);
    if (numFaces > 1) {
        std::ostringstream msg;
        msg << "Found multiple faces corresponding to vertices";
        for (PetscInt i = 0; i < numFaceVertices; ++i) {
            msg << " " << faceVertices[i];
        } // for
        msg << " in cell " << cell << ".";
        err = DMPlexRestoreJoin(dmMesh, numFaceVertices, faceVertices, &numFaces, &faces);PYLITH_CHECK_ERROR(err);
        throw std::runtime_error(msg.str());
    } else if (0 == numFaces) {
        std::ostringstream msg;
        msg << "Could not find face corresponding to vertices";
        for (PetscInt i = 0; i < numFaceVertices; ++i) {
            msg << " " << faceVertices[i];
        } // for
        msg << " in cell " << cell << ".";
        err = DMPlexRestoreJoin(dmMesh, numFaceVertices, faceVertices, &numFaces, &faces);PYLITH_CHECK_ERROR(err);
        throw std::runtime_error(msg.str());
    } // if
    *face = faces[0];
    err = DMPlexRestoreJoin(dmMesh, numFaceVertices, faceVertices, &numFaces, &faces);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // faceFromCellVertices


// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_MeshBuilder::cellVerticesFromFace(PetscInt* cell,
                                                   pylith::int_array* faceVertices,
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

    pylith::topology::Stratum verticesStratum(dmMesh, pylith::topology::Stratum::DEPTH, 0);
    const PetscInt vBegin = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    PetscInt closureSize = 0;
    PetscInt* closure = PETSC_NULLPTR;
    err = DMPlexGetTransitiveClosure(dmMesh, face, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    int_array buffer(closureSize);
    size_t numFaceVertices = 0;
    for (PetscInt iClosure = 0; iClosure < closureSize*2; iClosure += 2) {
        if ((closure[iClosure] >= vBegin) && (closure[iClosure] < vEnd)) {
            buffer[numFaceVertices++] = closure[iClosure];
        } // if
    } // for
    err = DMPlexRestoreTransitiveClosure(dmMesh, face, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    *faceVertices = int_array(&buffer[0], numFaceVertices);

    PYLITH_METHOD_END;
} // cellVerticesFromFace


// End of file
