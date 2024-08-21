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

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#include <algorithm> // USES std::sort, std::find
#include <map> // USES std::map

// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Create subdomain mesh using label.
pylith::topology::Mesh*
pylith::topology::MeshOps::createSubdomainMesh(const pylith::topology::Mesh& mesh,
                                               const char* labelName,
                                               const int labelValue,
                                               const char* descriptiveLabel) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createSubdomainMesh);

    assert(labelName);

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscErrorCode err = 0;

    PetscBool hasLabel = PETSC_FALSE;
    err = DMHasLabel(dmDomain, labelName, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* :TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmDomain, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, labelValue, &hasLabelValue);PYLITH_CHECK_ERROR(err);
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    err = MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDM dmSubdomain = NULL;
    err = DMPlexFilter(dmDomain, dmLabel, labelValue, PETSC_FALSE, PETSC_FALSE, NULL, &dmSubdomain);PYLITH_CHECK_ERROR(err);

    PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
    err = DMPlexGetMaxSizes(dmSubdomain, &maxConeSizeLocal, NULL);PYLITH_CHECK_ERROR(err);
    err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmSubdomain));PYLITH_CHECK_ERROR(err);

    if (maxConeSize <= 0) {
        err = DMDestroy(&dmSubdomain);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Error while creating mesh of subdomain. Subdomain mesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Set lengthscale
    PylithScalar lengthScale;
    err = DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmSubdomain, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    pylith::topology::Mesh* submesh = new pylith::topology::Mesh(true);assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());
    submesh->setDM(dmSubdomain, descriptiveLabel);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createSubdomainMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createSubdomainMesh


// ---------------------------------------------------------------------------------------------------------------------
// Create lower dimension mesh using label.
pylith::topology::Mesh*
pylith::topology::MeshOps::createLowerDimMesh(const pylith::topology::Mesh& mesh,
                                              const char* labelName,
                                              const int labelValue) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createLowerDimMesh);
    assert(labelName);

    if (mesh.getDimension() < 1) {
        throw std::logic_error("INTERNAL ERROR in MeshOps::createLowerDimMesh()\n"
                               "Cannot create submesh for mesh with dimension < 1.");
    } // if

    PetscErrorCode err = PETSC_SUCCESS;

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscBool hasLabel = PETSC_FALSE;
    err = DMHasLabel(dmDomain, labelName, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmDomain, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, labelValue, &hasLabelValue);PYLITH_CHECK_ERROR(err);
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    err = MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDM dmSubmesh = NULL;
    err = DMPlexCreateSubmesh(dmDomain, dmLabel, labelValue, PETSC_FALSE, &dmSubmesh);PYLITH_CHECK_ERROR(err);

    PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
    err = DMPlexGetMaxSizes(dmSubmesh, &maxConeSizeLocal, NULL);PYLITH_CHECK_ERROR(err);
    err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmSubmesh));PYLITH_CHECK_ERROR(err);

    if (maxConeSize <= 0) {
        err = DMDestroy(&dmSubmesh);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Error while creating lower dimension mesh. Submesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Set lengthscale
    PylithScalar lengthScale;
    err = DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmSubmesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    pylith::topology::Mesh* submesh = new pylith::topology::Mesh(true);assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());

    std::string meshLabel = "subdomain_" + std::string(labelName);
    submesh->setDM(dmSubmesh, meshLabel.c_str());

    // Check topology
    MeshOps::checkTopology(*submesh);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createLowerDimMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createLowerDimMesh


// ---------------------------------------------------------------------------------------------------------------------
// Create 0-dimension mesh from points.
pylith::topology::Mesh*
pylith::topology::MeshOps::createFromPoints(const PylithReal* points,
                                            const size_t numPoints,
                                            const spatialdata::geocoords::CoordSys* cs,
                                            const PylithReal lengthScale,
                                            MPI_Comm comm) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createFromPoints);
    assert(cs);

    PetscErrorCode err;

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

    err = DMPlexCreate(comm, &dmPoints);PYLITH_CHECK_ERROR(err);
    err = DMSetDimension(dmPoints, 0);PYLITH_CHECK_ERROR(err);
    err = DMSetCoordinateDim(dmPoints, spaceDim);PYLITH_CHECK_ERROR(err);
    if (numPoints > 0) {
        err = DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &dmConeSizes[0], &dmCones[0],
                                  &dmConeOrientations[0], points);PYLITH_CHECK_ERROR(err);
    } else {
        PetscInt empty[1];
        empty[0] = 0;
        err = DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &empty[0], &empty[0],
                                  &empty[0], points);PYLITH_CHECK_ERROR(err);
    } // if/else

    PetscSF sf = NULL;
    err = DMGetPointSF(dmPoints, &sf);PYLITH_CHECK_ERROR(err);
    err = PetscSFSetGraph(sf, numPoints, 0, NULL, PETSC_COPY_VALUES, NULL, PETSC_COPY_VALUES);

    mesh->setDM(dmPoints, "points");

    mesh->setCoordSys(cs);

    err = DMPlexSetScale(mesh->getDM(), PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createFromPoints);
    PYLITH_METHOD_RETURN(mesh);
} // createFromPoints


// ---------------------------------------------------------------------------------------------------------------------
// Nondimensionalize the finite-element mesh.
void
pylith::topology::MeshOps::nondimensionalize(Mesh* const mesh,
                                             const spatialdata::units::Nondimensional& normalizer) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::nondimensionalize);

    assert(mesh);

    PetscVec coordVec = NULL;
    const PylithScalar lengthScale = normalizer.getLengthScale();
    PetscErrorCode err = 0;

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);assert(coordVec);
    err = VecScale(coordVec, 1.0/lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dmMesh, NULL, "-pylith_nondim_dm_view");PYLITH_CHECK_ERROR(err);

    const PetscInt dim = mesh->getDimension();
    if (dim < 1) {
        PYLITH_METHOD_END;
    } // if
    PylithReal coordMin[3];
    PylithReal coordMax[3];
    err = DMGetBoundingBox(dmMesh, coordMin, coordMax);
    PylithReal volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= coordMax[i] - coordMin[i];
    } // for
    assert(dim > 0);
    const PylithReal avgCellDim = pow(volume / MeshOps::getNumCells(*mesh), 1.0/dim);
    const PylithReal avgDimTolerance = 0.02;
    if (avgCellDim < avgDimTolerance) {
        std::ostringstream msg;
        msg << "Nondimensional average cell dimension (" << avgCellDim << ") is less than minimum tolerance ("
            << avgDimTolerance << "). This usually means the length scale (" << lengthScale << ") used in the "
            << "nondimensionalization needs to be smaller. Based on the average cell size, a value of about "
            << pow(10, int(log10(avgCellDim*lengthScale))) << " should be appropriate.";
        throw std::runtime_error(msg.str());
    } // if/else

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::nondimensionalize);
    PYLITH_METHOD_END;
} // nondimensionalize


// ---------------------------------------------------------------------------------------------------------------------
// Strip out "ghost" cells hanging off mesh
PetscDM
pylith::topology::MeshOps::removeHangingCells(const PetscDM& dmMesh) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmClean = PETSC_NULLPTR;

    MPI_Comm comm = PetscObjectComm((PetscObject) dmMesh);
    pylith::topology::Stratum cells(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    DMPolytopeType cellType;
    err = DMPlexGetCellType(dmMesh, cells.begin(), &cellType);PYLITH_CHECK_ERROR(err);
    if (DMPolytopeTypeGetDim(cellType) < 0) {
        // Hanging cells have dim == -1

        // Create label over cells 1 dimension lower
        PetscDMLabel labelInclude = PETSC_NULLPTR;
        const PetscInt labelValue = 1;
        err = DMLabelCreate(comm, "no_hanging_cells", &labelInclude);PYLITH_CHECK_ERROR(err);
        pylith::topology::Stratum faces(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
        for (PetscInt face = faces.begin(); face < faces.end(); ++face) {
            err = DMLabelSetValue(labelInclude, face, labelValue);PYLITH_CHECK_ERROR(err);
        } // for

        err = DMPlexFilter(dmMesh, labelInclude, labelValue, PETSC_FALSE, PETSC_FALSE, PETSC_NULLPTR, &dmClean);PYLITH_CHECK_ERROR(err);
    } else {
        dmClean = dmMesh;
        err = PetscObjectReference((PetscObject) dmClean);
    } // if/else
    PetscReal lengthScale = 1.0;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmClean, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(dmClean);
}


// ---------------------------------------------------------------------------------------------------------------------
// Check topology of mesh.
void
pylith::topology::MeshOps::checkTopology(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopology);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    DMLabel subpointMap;
    PetscErrorCode ierr = DMPlexGetSubpointMap(dmMesh, &subpointMap);PYLITH_CHECK_ERROR(ierr);
    PetscInt cellHeight = subpointMap ? 1 : 0;

    PetscErrorCode err;
    err = DMViewFromOptions(dmMesh, NULL, "-pylith_checktopo_dm_view");PYLITH_CHECK_ERROR(err);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologyGeometry);
    err = DMPlexCheckGeometry(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of the mesh.");
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologyGeometry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySymmetry);
    err = DMPlexCheckSymmetry(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh associated with symmetry of adjacency information.");
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologySymmetry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySkeleton);
    err = DMPlexCheckSkeleton(dmMesh, cellHeight);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh cells.");
    err = DMPlexCheckOrphanVertices(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Mesh contains vertices not connected to cells.");
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


// ---------------------------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isSimplexMesh(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    bool isSimplex = false;

    PetscErrorCode err;
    const PetscDM dm = mesh.getDM();
    PetscInt vStart = 0, vEnd = 0;
    err = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    if (vStart != vEnd) { // Test for simplex only works if we have points.
        PetscInt closureSize = 0;
        PetscInt* closure = NULL;
        const int dim = mesh.getDimension();

        err = DMPlexGetTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        PetscInt numVertices = 0;
        for (PetscInt c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                ++numVertices;
            } // if
        } // for
        if (numVertices == dim+1) {
            isSimplex = PETSC_TRUE;
        } // if
        err = DMPlexRestoreTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // if

    // Communicate result of isSimplex to all processes.
    int intSimplexLocal = isSimplex ? 1 : 0;
    int intSimplexGlobal = 0;
    MPI_Allreduce(&intSimplexLocal, &intSimplexGlobal, 1, MPI_INT, MPI_LOR, mesh.getComm());
    isSimplex = intSimplexGlobal == 1;

    PYLITH_METHOD_RETURN(isSimplex);
} // isSimplexMesh


// ---------------------------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isCohesiveCell(const PetscDM dm,
                                          const PetscInt cell) {
    bool isCohesive = false;

    DMPolytopeType ct;
    PetscErrorCode err = DMPlexGetCellType(dm, cell, &ct);PYLITH_CHECK_ERROR(err);
    if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_TRI_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_QUAD_PRISM_TENSOR)) { isCohesive = true; }

    return isCohesive;
} // isCohesiveCell


// ----------------------------------------------------------------------
// Get number of vertices in mesh.
PylithInt
pylith::topology::MeshOps::getNumVertices(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PylithInt nvertices = 0;
    PylithInt begin = 0, end = 0;
    PetscErrorCode err = DMPlexGetDepthStratum(dmMesh, 0, &begin, &end);PYLITH_CHECK_ERROR(err);
    nvertices = end-begin;

    PYLITH_METHOD_RETURN(nvertices);
}


// ----------------------------------------------------------------------
// Get number of cells in mesh.
PylithInt
pylith::topology::MeshOps::getNumCells(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscInt ncells = 0;
    PylithInt begin = 0, end = 0;
    const int cellHeight = 0;
    PetscErrorCode err = DMPlexGetHeightStratum(dmMesh, cellHeight, &begin, &end);PYLITH_CHECK_ERROR(err);
    ncells = end-begin;

    PYLITH_METHOD_RETURN(ncells);
}


// ----------------------------------------------------------------------
// Get number of vertices in a cell.
PylithInt
pylith::topology::MeshOps::getNumCorners(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscInt numCorners = 0;
    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    PetscInt cStart, cEnd, vStart, vEnd, closureSize, *closure = NULL;
    PetscErrorCode err;
    const int cellHeight = 0;
    err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    if (cEnd > cStart) {
        err = DMPlexGetTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PetscInt c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {++numCorners;}
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(numCorners);
}


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialLabels(const pylith::topology::Mesh& mesh,
                                               pylith::int_array& labelValues) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkMaterialLabels);

    PetscErrorCode err;

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
    err = DMGetLabel(dmMesh, labelName, &materialsLabel);PYLITH_CHECK_ERROR(err);assert(materialsLabel);

    int *matBegin = &labelValues[0];
    int *matEnd = &labelValues[0] + labelValues.size();
    std::sort(matBegin, matEnd);

    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt matId;

        err = DMLabelGetValue(materialsLabel, c, &matId);PYLITH_CHECK_ERROR(err);
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
    err = MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
                        matCellCounts.size(), MPI_INT, MPI_SUM, mesh.getComm());PYLITH_CHECK_ERROR(err);
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
