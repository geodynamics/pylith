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

#include "pylith/meshio/MeshIOCubit.hh" // implementation of class methods

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/meshio/ExodusII.hh" // USES ExodusII

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "petsc.h" // USES MPI_Comm

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        namespace _MeshIOCubit {
            /** Read mesh vertices.
             *
             * @param[out] geometry Mesh geometry.
             * @param[inout] fileIn Cubit Exodus file.
             */
            void readVertices(pylith::meshio::MeshBuilder::Geometry* geometry,
                              pylith::meshio::ExodusII& fileIn);

            /** Write mesh vertices.
             *
             * @param[inout] fileOut Cubit Exodus file.
             * @param[in] geometry Mesh geometry.
             */
            void writeVertices(pylith::meshio::ExodusII& fileOut,
                               const pylith::meshio::MeshBuilder::Geometry& geometry);

            /** Read mesh cells.
             *
             * @param[out] topology Mesh topology.
             * @param[out] materialIds Material id for each cell.
             * @param[inout] fileIn Cubit Exodus file.
             */
            void readCells(pylith::meshio::MeshBuilder::Topology* topology,
                           int_array* materialIds,
                           pylith::meshio::ExodusII& fileIn);

            /** Write mesh cells.
             *
             * @param[inout] fileOut Cubit Exodus file.
             * @param[in] topology Mesh topology.
             */
            void writeCells(pylith::meshio::ExodusII& fileOut,
                            const pylith::meshio::MeshBuilder::Topology& topology,
                            const int_array& materialIds);

            /** Read a point group with vertices.
             *
             * @param[out] points Vertices in group.
             * @param[out] name Name of group.
             * @param[inout] fileIn Cubit Exodus file.
             */
            void readNodeSet(int_array* points,
                             std::string* name,
                             pylith::meshio::ExodusII& fileIn);

            /** Write a point group with vertices.
             *
             * @param[inout] fileOut Output stream.
             * @param[in] points Vertices in group.
             * @param[in] name Name of group.
             */
            void writeNodeSet(std::ostream& fileOut,
                              const int_array& points,
                              const char* name);

            /** Read a point group with faces.
             *
             * @param[out] faceValues Array of cell+vertices for each face.
             * @param[out] name Name of group.
             * @param[inout] parser Input parser.
             * @param[in] faceShape Shape of face.
             * @param[in] useIndexZero True if using zero-based indexing.
             */
            void readSideSet(int_array* faceValues,
                             std::string* name,
                             pylith::meshio::ExodusII& fileIn);

            /** Write a point group with faces
             *
             * @param[inout] fileOut Output stream.
             * @param[in] faceValues Array of cell+vertices for each face.
             * @param[in] faceShape Shape of face.
             * @param[in] name Name of group.
             */
            void writeSideSet(std::ostream& fileOut,
                              const int_array& faceValues,
                              const char* name);

            /** Reorder vertices in cells to match PyLith/PETSc conventions.
             *
             * @param[inout] topology Mesh topology.
             */
            void orientCells(pylith::meshio::MeshBuilder::Topology* topology);

        } // _MeshIOCubit
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOCubit::MeshIOCubit(void) :
    _filename("") {
    PyreComponent::setName("meshiocubit");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOCubit::~MeshIOCubit(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOCubit::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set filename for CUBIT file.
void
pylith::meshio::MeshIOCubit::setFilename(const char* name) {
    _filename = name;
} // setFilename


// ------------------------------------------------------------------------------------------------
// Get filename of CUBIT file.
const char*
pylith::meshio::MeshIOCubit::getFilename(void) const {
    return _filename.c_str();
} // getFilename


// ------------------------------------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOCubit::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    assert(_mesh);

    const int commRank = _mesh->getCommRank();

    pylith::meshio::MeshBuilder::Topology topology;
    pylith::meshio::MeshBuilder::Geometry geometry;
    int_array materialIds;

    if (0 == commRank) {
        try {
            pylith::meshio::ExodusII exoFile(_filename.c_str());

            topology.dimension = exoFile.getDim("num_dim");

            _MeshIOCubit::readVertices(&geometry, exoFile);
            PYLITH_COMPONENT_INFO_ROOT("Read " << geometry.numVertices << " vertices.");

            _MeshIOCubit::readCells(&topology, &materialIds, exoFile);
            topology.cellShape = pylith::meshio::MeshBuilder::cellShapeFromCorners(topology.dimension, topology.numCorners);
            PYLITH_COMPONENT_INFO_ROOT("Read " << topology.numCells << " cells.");

            _MeshIOCubit::orientCells(&topology);
            pylith::meshio::MeshBuilder::buildMesh(_mesh, topology, geometry);
            pylith::meshio::MeshBuilder::setMaterials(_mesh, materialIds);

            _readNodeSets(exoFile);
            _readSideSets(exoFile);
        } catch (std::exception& err) {
            std::ostringstream msg;
            msg << "Error while reading Cubit Exodus file '" << _filename << "'.\n"
                << err.what();
            throw std::runtime_error(msg.str());
        } catch (...) {
            std::ostringstream msg;
            msg << "Unknown error while reading Cubit Exodus file '" << _filename << "'.";
            throw std::runtime_error(msg.str());
        } // try/catch
    } else {
        pylith::meshio::MeshBuilder::buildMesh(_mesh, topology, geometry);
        pylith::meshio::MeshBuilder::setMaterials(_mesh, materialIds);
    } // if/else

    PYLITH_METHOD_END;
} // read


// ------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOCubit::_write(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_write()");

    throw std::logic_error("MeshIOCubit::_write() not implemented.");

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshio::_MeshIOCubit::readVertices(pylith::meshio::MeshBuilder::Geometry* geometry,
                                           pylith::meshio::ExodusII& exoFile) {
    PYLITH_METHOD_BEGIN;

    assert(geometry);

    geometry->spaceDim = exoFile.getDim("num_dim");
    geometry->numVertices = exoFile.getDim("num_nodes");

    if (exoFile.hasVar("coord", NULL)) {
        const int ndims = 2;
        int dims[2];
        dims[0] = geometry->spaceDim;
        dims[1] = geometry->numVertices;
        scalar_array buffer(geometry->numVertices * geometry->spaceDim);
        exoFile.getVar(&buffer[0], dims, ndims, "coord");

        geometry->vertices.resize(geometry->numVertices * geometry->spaceDim);
        for (size_t iVertex = 0; iVertex < geometry->numVertices; ++iVertex) {
            for (size_t iDim = 0; iDim < geometry->spaceDim; ++iDim) {
                geometry->vertices[iVertex*(geometry->spaceDim)+iDim] =
                    buffer[iDim*geometry->numVertices+iVertex];
            } // for
        } // for

    } else {
        const char* coordNames[3] = { "coordx", "coordy", "coordz" };

        geometry->vertices.resize(geometry->numVertices * geometry->spaceDim);
        scalar_array buffer(geometry->numVertices);

        const int ndims = 1;
        int dims[1];
        dims[0] = geometry->numVertices;

        for (size_t i = 0; i < geometry->spaceDim; ++i) {
            exoFile.getVar(&buffer[0], dims, ndims, coordNames[i]);

            for (size_t iVertex = 0; iVertex < geometry->numVertices; ++iVertex) {
                geometry->vertices[iVertex*(geometry->spaceDim)+i] = buffer[iVertex];
            } // for
        } // for
    } // else

    PYLITH_METHOD_END;
} // readVertices


// ------------------------------------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshio::_MeshIOCubit::readCells(pylith::meshio::MeshBuilder::Topology* topology,
                                        int_array* materialIds,
                                        pylith::meshio::ExodusII& exoFile) {
    PYLITH_METHOD_BEGIN;

    assert(topology);
    assert(materialIds);

    topology->numCells = exoFile.getDim("num_elem");
    const int numMaterials = exoFile.getDim("num_el_blk");

    int_array blockIds(numMaterials);
    int ndims = 1;
    int dims[2];
    dims[0] = numMaterials;
    dims[1] = 0;
    exoFile.getVar(&blockIds[0], dims, ndims, "eb_prop1");

    materialIds->resize(topology->numCells);
    topology->numCorners = 0;
    for (int iMaterial = 0, index = 0; iMaterial < numMaterials; ++iMaterial) {
        std::ostringstream varName;
        varName << "num_nod_per_el" << iMaterial+1;
        if (0 == topology->numCorners) {
            topology->numCorners = exoFile.getDim(varName.str().c_str());
            const int size = (topology->numCells) * (topology->numCorners);
            topology->cells.resize(size);
        } else if (size_t(exoFile.getDim(varName.str().c_str())) != topology->numCorners) {
            std::ostringstream msg;
            msg << "All materials must have the same number of vertices per cell.\n"
                << "Expected " << topology->numCorners << " vertices per cell, but block "
                << blockIds[iMaterial] << " has "
                << exoFile.getDim(varName.str().c_str())
                << " vertices.";
            throw std::runtime_error(msg.str());
        } // if

        varName.str("");
        varName << "num_el_in_blk" << iMaterial+1;
        const int blockSize = exoFile.getDim(varName.str().c_str());

        varName.str("");
        varName << "connect" << iMaterial+1;
        ndims = 2;
        dims[0] = blockSize;
        dims[1] = topology->numCorners;
        exoFile.getVar(&topology->cells[index* topology->numCorners], dims, ndims, varName.str().c_str());

        for (int i = 0; i < blockSize; ++i) {
            (*materialIds)[index+i] = blockIds[iMaterial];
        } // for

        index += blockSize;
    } // for

    topology->cells -= 1; // use zero index

    PYLITH_METHOD_END;
} // _readCells


// ------------------------------------------------------------------------------------------------
// Read vertex groups.
void
pylith::meshio::MeshIOCubit::_readNodeSets(ExodusII& exoFile) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_readNodeSets(exoFile="<<typeid(exoFile).name()<<")");

    if (!exoFile.hasDim("num_node_sets", NULL)) {
        PYLITH_COMPONENT_INFO_ROOT("No nodesets found.");
        PYLITH_METHOD_END;
    } // if
    const int numGroups = exoFile.getDim("num_node_sets");
    PYLITH_COMPONENT_INFO_ROOT("Found " << numGroups << " node sets.");
    if (!numGroups) {
        PYLITH_METHOD_END;
    } // if

    int_array ids(numGroups);
    int ndims = 1;
    int dims[2];
    dims[0] = numGroups;
    dims[1] = 0;
    exoFile.getVar(&ids[0], dims, ndims, "ns_prop1");

    string_vector groupNames(numGroups);
    exoFile.getVar(&groupNames, numGroups, "ns_names");

    for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
        std::ostringstream varName;
        varName << "num_nod_ns" << iGroup+1;
        const size_t nodesetSize = exoFile.getDim(varName.str().c_str());
        int_array points(nodesetSize);

        varName.str("");
        varName << "node_ns" << iGroup+1;
        ndims = 1;
        dims[0] = nodesetSize;

        PYLITH_COMPONENT_INFO_ROOT("Reading node set '" << groupNames[iGroup] << "' with id " << ids[iGroup] << " containing " << nodesetSize << " nodes.");
        exoFile.getVar(&points[0], dims, ndims, varName.str().c_str());

        std::sort(&points[0], &points[0]+nodesetSize);
        points -= 1; // use zero index

        pylith::meshio::MeshBuilder::setVertexGroup(_mesh, groupNames[iGroup].c_str(), points);
    } // for

    PYLITH_METHOD_END;
} // _readNodeSets


// ------------------------------------------------------------------------------------------------
// Read face groups.
void
pylith::meshio::MeshIOCubit::_readSideSets(ExodusII& exoFile) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_readSideSets(exoFile="<<typeid(exoFile).name()<<")");

    if (!exoFile.hasDim("num_side_sets", NULL)) {
        PYLITH_COMPONENT_INFO_ROOT("No sidesets found.");
        PYLITH_METHOD_END;
    } // if
    const int numGroups = exoFile.getDim("num_side_sets");
    PYLITH_COMPONENT_INFO_ROOT("Found " << numGroups << " side sets.");
    if (!numGroups) {
        PYLITH_METHOD_END;
    } // if

    // Need cell shape to interpret sides
    const char* varName = "connect1"; // :KLUDGE: Assume uniform cell type
    const char* attrName = "elem_type";
    std::string shapeName;
    exoFile.getAttr(&shapeName, varName, attrName);
    int sideOffset = -1; // 1-based to 0-based index
    if ((shapeName == std::string("SHELL")) || (shapeName == std::string("SHELL4"))) {
        sideOffset -= 2; // 6 sides and start at 3 instead of 1
    } // if

    int_array ids(numGroups);
    int ndims = 1;
    int dims[2];
    dims[0] = numGroups;
    dims[1] = 0;
    exoFile.getVar(&ids[0], dims, ndims, "ss_prop1");

    string_vector groupNames(numGroups);
    exoFile.getVar(&groupNames, numGroups, "ss_names");

    for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
        std::ostringstream varName;
        varName << "num_side_ss" << iGroup+1;
        const size_t sideSetSize = exoFile.getDim(varName.str().c_str());
        int_array ioBuffer(sideSetSize);
        int_array points(2*sideSetSize);

        PYLITH_COMPONENT_INFO_ROOT("Reading side set '" << groupNames[iGroup] << "' with id " << ids[iGroup] << " containing " << sideSetSize << " faces.");

        // Read cells
        varName.str("");
        varName << "elem_ss" << iGroup+1;
        ndims = 1;
        dims[0] = sideSetSize;
        exoFile.getVar(&ioBuffer[0], dims, ndims, varName.str().c_str());
        ioBuffer -= 1; // use zero index
        for (size_t i = 0; i < sideSetSize; ++i) {
            points[2*i+0] = ioBuffer[i];
        } // for

        // Read cell sides
        varName.str("");
        varName << "side_ss" << iGroup+1;
        ndims = 1;
        dims[0] = sideSetSize;
        exoFile.getVar(&ioBuffer[0], dims, ndims, varName.str().c_str());
        ioBuffer += sideOffset;
        for (size_t i = 0; i < sideSetSize; ++i) {
            points[2*i+1] = ioBuffer[i];
        } // for

        pylith::meshio::MeshBuilder::setFaceGroupFromCellSide(_mesh, groupNames[iGroup].c_str(), points);
    } // for

    PYLITH_METHOD_END;
} // _readSideSets


// ------------------------------------------------------------------------------------------------
// Reorder vertices in cells to match PyLith conventions.
void
pylith::meshio::_MeshIOCubit::orientCells(pylith::meshio::MeshBuilder::Topology* topology) {
    PYLITH_METHOD_BEGIN;

    assert(topology);
    const size_t numCorners = topology->numCorners;

    if ((2 == topology->dimension) && (6 == numCorners)) { // TRI6
        // CUBIT
        // corners,
        // bottom edges, middle edges, top edges

        // PyLith
        // bottom edge, right edge, left edge, corners

        // Permutation: 3, 4, 5, 0, 1, 2
        int tmp = 0;
        for (size_t iCell = 0; iCell < topology->numCells; ++iCell) {
            const size_t ii = iCell*numCorners;
            tmp = (topology->cells)[ii+0];
            (topology->cells)[ii+0] = (topology->cells)[ii+3];
            (topology->cells)[ii+3] = tmp;

            tmp = (topology->cells)[ii+1];
            (topology->cells)[ii+1] = (topology->cells)[ii+4];
            (topology->cells)[ii+4] = tmp;

            tmp = (topology->cells)[ii+2];
            (topology->cells)[ii+2] = (topology->cells)[ii+5];
            (topology->cells)[ii+5] = tmp;
        } // for

    } else if ((3 == topology->dimension) && (27 == numCorners)) { // HEX27
        // CUBIT
        // corners,
        // bottom edges, middle edges, top edges
        // interior
        // bottom/top, left/right, front/back

        // PyLith
        // corners,
        // bottom edges, top edges, middle edges
        // left/right, front/back, bottom/top
        // interior
        int tmp = 0;
        for (size_t iCell = 0; iCell < topology->numCells; ++iCell) {
            const int i12 = iCell*numCorners+12;
            const int i13 = iCell*numCorners+13;
            const int i14 = iCell*numCorners+14;
            const int i15 = iCell*numCorners+15;
            const int i16 = iCell*numCorners+16;
            const int i17 = iCell*numCorners+17;
            const int i18 = iCell*numCorners+18;
            const int i19 = iCell*numCorners+19;
            const int i20 = iCell*numCorners+20;
            const int i21 = iCell*numCorners+21;
            const int i22 = iCell*numCorners+22;
            const int i23 = iCell*numCorners+23;
            const int i24 = iCell*numCorners+24;
            const int i25 = iCell*numCorners+25;
            const int i26 = iCell*numCorners+26;

            tmp = (topology->cells)[i12];
            (topology->cells)[i12] = (topology->cells)[i16];
            (topology->cells)[i16] = tmp;

            tmp = (topology->cells)[i13];
            (topology->cells)[i13] = (topology->cells)[i17];
            (topology->cells)[i17] = tmp;

            tmp = (topology->cells)[i14];
            (topology->cells)[i14] = (topology->cells)[i18];
            (topology->cells)[i18] = tmp;

            tmp = (topology->cells)[i15];
            (topology->cells)[i15] = (topology->cells)[i19];
            (topology->cells)[i19] = tmp;

            tmp = (topology->cells)[i20];
            (topology->cells)[i20] = (topology->cells)[i23];
            (topology->cells)[i23] = (topology->cells)[i26];
            (topology->cells)[i26] = tmp;

            tmp = (topology->cells)[i21];
            (topology->cells)[i21] = (topology->cells)[i24];
            (topology->cells)[i24] = tmp;

            tmp = (topology->cells)[i22];
            (topology->cells)[i22] = (topology->cells)[i25];
            (topology->cells)[i25] = tmp;
        } // for
    } // if/else

    PYLITH_METHOD_END;
} // _orientCells


// End of file
