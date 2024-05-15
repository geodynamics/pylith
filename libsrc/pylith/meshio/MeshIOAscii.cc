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

#include "pylith/meshio/MeshIOAscii.hh" // implementation of class methods

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "spatialdata/utils/LineParser.hh" // USES LineParser

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iomanip> // USES setw(), setiosflags(), resetiosflags()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        namespace _MeshIOAscii {
            /** Read mesh vertices.
             *
             * @param[out] geometry Mesh geometry.
             * @param[inout] parser Input parser.
             */
            void readVertices(pylith::meshio::MeshBuilder::Geometry* geometry,
                              spatialdata::utils::LineParser& parser);

            /** Write mesh vertices.
             *
             * @param[inout[ fileout Output stream
             * @param[in] geometry Mesh geometry.
             */
            void writeVertices(std::ostream& fileout,
                               const pylith::meshio::MeshBuilder::Geometry& geometry);

            /** Read mesh cells.
             *
             * @param[out] topology Mesh topology.
             * @param[out] materialIds Material id for each cell.
             * @param[inout] parser Input parser.
             * @param[in] useIndexZero True if using zero-based indexing.
             */
            void readCells(pylith::meshio::MeshBuilder::Topology* topology,
                           int_array* materialIds,
                           spatialdata::utils::LineParser& parser,
                           const bool useIndexZero);

            /** Write mesh cells.
             *
             * @param fileout Output stream
             * @param[in] topology Mesh topology.
             */
            void writeCells(std::ostream& fileout,
                            const pylith::meshio::MeshBuilder::Topology& topology,
                            const int_array& materialIds);

            /** Read a point group with vertices.
             *
             * @param[out] points Vertices in group.
             * @param[out] name Name of group.
             * @param[inout] parser Input parser.
             * @param[in] useIndexZero True if using zero-based indexing.
             */
            void readVertexGroup(int_array* points,
                                 std::string* name,
                                 spatialdata::utils::LineParser& parser,
                                 const bool useIndexZero);

            /** Write a point group with vertices.
             *
             * @param[inout] fileout Output stream.
             * @param[in] points Vertices in group.
             * @param[in] name Name of group.
             */
            void writeVertexGroup(std::ostream& fileout,
                                  const int_array& points,
                                  const char* name);

            /** Read a point group with facs.
             *
             * @param[out] faceValues Array of cell+vertices for each face.
             * @param[out] name Name of group.
             * @param[inout] parser Input parser.
             * @param[in] faceShape Shape of face.
             * @param[in] useIndexZero True if using zero-based indexing.
             */
            void readFaceGroup(int_array* faceValues,
                               std::string* name,
                               spatialdata::utils::LineParser& parser,
                               pylith::meshio::MeshBuilder::shape_t faceShape,
                               const bool useIndexZero);

            /** Write a point group with faces
             *
             * @param[inout] fileout Output stream.
             * @param[in] faceValues Array of cell+vertices for each face.
             * @param[in] faceShape Shape of face.
             * @param[in] name Name of group.
             */
            void writeFaceGroup(std::ostream& fileout,
                                const int_array& faceValues,
                                const pylith::meshio::MeshBuilder::shape_t faceShape,
                                const char* name);

        } // _MeshIOAscii
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOAscii::MeshIOAscii(void) :
    _filename(""),
    _useIndexZero(true) {
    PyreComponent::setName("meshioascii");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOAscii::~MeshIOAscii(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOAscii::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set filename for ASCII file.
void
pylith::meshio::MeshIOAscii::setFilename(const char* name) {
    _filename = name;
}


// ------------------------------------------------------------------------------------------------
// Get filename of ASCII file.
const char*
pylith::meshio::MeshIOAscii::getFilename(void) const {
    return _filename.c_str();
}


// ------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOAscii::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");

    const int commRank = _mesh->getCommRank();
    pylith::meshio::MeshBuilder::Topology topology;
    pylith::meshio::MeshBuilder::Geometry geometry;
    int_array materialIds;

    if (0 == commRank) {
        std::ifstream filein(_filename.c_str());
        if (!filein.is_open() || !filein.good()) {
            std::ostringstream msg;
            msg << "Could not open mesh file '" << _filename
                << "' for reading.\n";
            throw std::runtime_error(msg.str());
        } // if

        spatialdata::utils::LineParser parser(filein, "//");
        parser.eatwhitespace(true);

        std::string token;
        std::istringstream buffer;
        const int maxIgnore = 1024;

        buffer.str(parser.next());
        buffer >> token;
        if (strcasecmp(token.c_str(), "mesh")) {
            std::ostringstream msg;
            msg << "Expected 'mesh' token but encountered '" << token << "'\n";
            throw std::runtime_error(msg.str());
        } // if

        bool readDim = false;
        bool readCells = false;
        bool readVertices = false;
        bool builtMesh = false;

        try {
            buffer.str(parser.next());
            buffer.clear();
            buffer >> token;
            while (buffer.good() && token != "}") {
                if (0 == strcasecmp(token.c_str(), "dimension")) {
                    buffer.ignore(maxIgnore, '=');
                    buffer >> topology.dimension;
                    readDim = true;
                } else if (0 == strcasecmp(token.c_str(), "use-index-zero")) {
                    buffer.ignore(maxIgnore, '=');
                    std::string flag = "";
                    buffer >> flag;
                    if (0 == strcasecmp(flag.c_str(), "true")) {
                        _useIndexZero = true;
                    } else {
                        _useIndexZero = false;
                    }
                } else if (0 == strcasecmp(token.c_str(), "vertices")) {
                    _MeshIOAscii::readVertices(&geometry, parser);
                    readVertices = true;
                } else if (0 == strcasecmp(token.c_str(), "cells")) {
                    _MeshIOAscii::readCells(&topology, &materialIds, parser, _useIndexZero);
                    readCells = true;
                } else if (0 == strcasecmp(token.c_str(), "vertex-group")) {
                    if (!builtMesh) {
                        throw std::runtime_error("Both 'vertices' and 'cells' must "
                                                 "precede any groups in the mesh file.");
                    } // if

                    std::string name;
                    int_array points;
                    _MeshIOAscii::readVertexGroup(&points, &name, parser, _useIndexZero);
                    pylith::meshio::MeshBuilder::setVertexGroup(_mesh, name.c_str(), points);
                } else if (0 == strcasecmp(token.c_str(), "face-group")) {
                    if (!builtMesh) {
                        throw std::runtime_error("Both 'vertices' and 'cells' must "
                                                 "precede any groups in the mesh file.");
                    } // if

                    std::string name;
                    int_array faceValues;
                    pylith::meshio::MeshBuilder::shape_t faceShape = pylith::meshio::MeshBuilder::faceShapeFromCellShape(topology.cellShape);
                    _MeshIOAscii::readFaceGroup(&faceValues, &name, parser, faceShape, _useIndexZero);
                    pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(_mesh, name.c_str(), faceValues, faceShape);
                } else {
                    std::ostringstream msg;
                    msg << "Could not parse '" << token << "' into a mesh setting.";
                    throw std::runtime_error(msg.str());
                } // else

                if (readDim && readCells && readVertices && !builtMesh) {
                    // Can now build mesh
                    MeshBuilder::buildMesh(_mesh, topology, geometry);
                    pylith::meshio::MeshBuilder::setMaterials(_mesh, materialIds);
                    builtMesh = true;
                } // if

                buffer.str(parser.next());
                buffer.clear();
                buffer >> token;
            } // while
            if (token != "}") {
                throw std::runtime_error("I/O error occurred while parsing mesh tokens.");
            }
        } catch (const std::exception& err) {
            std::ostringstream msg;
            msg << "Error occurred while reading PyLith mesh ASCII file '"
                << _filename << "'.\n"
                << err.what();
            throw std::runtime_error(msg.str());
        } catch (...) {
            std::ostringstream msg;
            msg << "Unknown I/O error while reading PyLith mesh ASCII file '"
                << _filename << "'.\n";
            throw std::runtime_error(msg.str());
        } // catch
        filein.close();
    } else {
        pylith::meshio::MeshBuilder::buildMesh(_mesh, topology, geometry);
        pylith::meshio::MeshBuilder::setMaterials(_mesh, materialIds);
    } // if/else

    PYLITH_METHOD_END;
} // read


// ------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOAscii::_write(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_write()");

    std::ofstream fileout(_filename.c_str());
    if (!fileout.is_open() || !fileout.good()) {
        std::ostringstream msg;
        msg << "Could not open mesh file '" << _filename
            << "' for writing.\n";
        throw std::runtime_error(msg.str());
    } // if

    fileout
        << "mesh = {\n"
        << "  dimension = " << _mesh->getDimension() << "\n"
        << "  use-index-zero = " << (_useIndexZero ? "true" : "false") << "\n";

    pylith::meshio::MeshBuilder::Geometry geometry;
    pylith::meshio::MeshBuilder::getVertices(&geometry, *_mesh);
    _MeshIOAscii::writeVertices(fileout, geometry);

    pylith::meshio::MeshBuilder::Topology topology;
    pylith::meshio::MeshBuilder::getCells(&topology, *_mesh);
    int_array materialIds;
    pylith::meshio::MeshBuilder::getMaterials(&materialIds, *_mesh);
    _MeshIOAscii::writeCells(fileout, topology, materialIds);

    string_vector groupNames;

    pylith::meshio::MeshBuilder::getVertexGroupNames(&groupNames, *_mesh);
    size_t numGroups = groupNames.size();
    for (int i = 0; i < numGroups; ++i) {
        int_array points;
        pylith::meshio::MeshBuilder::getVertexGroup(&points, *_mesh, groupNames[i].c_str());
        _MeshIOAscii::writeVertexGroup(fileout, points, groupNames[i].c_str());
    } // for

    pylith::meshio::MeshBuilder::getFaceGroupNames(&groupNames, *_mesh);
    numGroups = groupNames.size();
    for (int i = 0; i < numGroups; ++i) {
        int_array faceValues;
        pylith::meshio::MeshBuilder::shape_t faceShape = pylith::meshio::MeshBuilder::faceShapeFromCellShape(topology.cellShape);
        pylith::meshio::MeshBuilder::getFaceGroup(&faceValues, *_mesh, groupNames[i].c_str());
        _MeshIOAscii::writeFaceGroup(fileout, faceValues, faceShape, groupNames[i].c_str());
    } // for

    fileout << "}\n";
    fileout.close();

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshio::_MeshIOAscii::readVertices(pylith::meshio::MeshBuilder::Geometry* geometry,
                                           spatialdata::utils::LineParser& parser) {
    PYLITH_METHOD_BEGIN;

    assert(geometry);

    std::string token;
    std::istringstream buffer;
    const int maxIgnore = 1024;
    buffer.str(parser.next());
    buffer.clear();
    buffer >> token;
    while (buffer.good() && token != "}") {
        if (0 == strcasecmp(token.c_str(), "dimension")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> geometry->spaceDim;
        } else if (0 == strcasecmp(token.c_str(), "count")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> geometry->numVertices;
        } else if (0 == strcasecmp(token.c_str(), "coordinates")) {
            const int size = (geometry->numVertices) * (geometry->spaceDim);
            if (0 == size) {
                const char* msg =
                    "Tokens 'dimension' and 'count' must precede 'coordinates'.";
                throw std::runtime_error(msg);
            } // if
            geometry->vertices.resize(size);
            int label;
            for (int iVertex = 0, i = 0; iVertex < geometry->numVertices; ++iVertex) {
                buffer.str(parser.next());
                buffer.clear();
                buffer >> label;
                for (int iDim = 0; iDim < geometry->spaceDim; ++iDim) {
                    buffer >> geometry->vertices[i++];
                }
            } // for
            parser.ignore('}');
        } else {
            std::ostringstream msg;
            msg << "Could not parse '" << token << "' into a vertices setting.";
            throw std::runtime_error(msg.str());
        } // else
        buffer.str(parser.next());
        buffer.clear();
        buffer >> token;
    } // while
    if (token != "}") {
        throw std::runtime_error("I/O error while parsing vertices.");
    }

    PYLITH_METHOD_END;
} // readVertices


// ------------------------------------------------------------------------------------------------
// Write mesh vertices.
void
pylith::meshio::_MeshIOAscii::writeVertices(std::ostream& fileout,
                                            const pylith::meshio::MeshBuilder::Geometry& geometry) {
    PYLITH_METHOD_BEGIN;

    fileout
        << "  vertices = {\n"
        << "    dimension = " << geometry.spaceDim << "\n"
        << "    count = " << geometry.numVertices << "\n"
        << "    coordinates = {\n"
        << std::resetiosflags(std::ios::fixed)
        << std::setiosflags(std::ios::scientific)
        << std::setprecision(6);
    for (int iVertex = 0, i = 0; iVertex < geometry.numVertices; ++iVertex) {
        fileout << "      ";
        fileout << std::setw(8) << iVertex;
        for (int iDim = 0; iDim < geometry.spaceDim; ++iDim) {
            fileout << std::setw(18) << geometry.vertices[i++];
        }
        fileout << "\n";
    } // for
    fileout
        << "    }\n"
        << "  }\n";

    PYLITH_METHOD_END;
} // writeVertices


// ------------------------------------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshio::_MeshIOAscii::readCells(pylith::meshio::MeshBuilder::Topology* topology,
                                        int_array* materialIds,
                                        spatialdata::utils::LineParser& parser,
                                        const bool useIndexZero) {
    PYLITH_METHOD_BEGIN;

    assert(topology);
    assert(materialIds);

    std::string token;
    std::istringstream buffer;
    const int maxIgnore = 1024;
    buffer.str(parser.next());
    buffer.clear();
    buffer >> token;
    while (buffer.good() && token != "}") {
        if (0 == strcasecmp(token.c_str(), "num-corners")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> topology->numCorners;
        } else if (0 == strcasecmp(token.c_str(), "count")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> topology->numCells;
        } else if (0 == strcasecmp(token.c_str(), "simplices")) {
            const int size = (topology->numCells) * (topology->numCorners);
            if (0 == size) {
                const char* msg =
                    "Tokens 'num-corners' and 'count' must precede 'cells'.";
                throw std::runtime_error(msg);
            } // if
            topology->cells.resize(size);
            int label;
            for (int iCell = 0, i = 0; iCell < topology->numCells; ++iCell) {
                buffer.str(parser.next());
                buffer.clear();
                buffer >> label;
                for (int iCorner = 0; iCorner < topology->numCorners; ++iCorner) {
                    buffer >> topology->cells[i++];
                }
            } // for
            if (!useIndexZero) {
                // if files begins with index 1, then decrement to index 0
                // for compatibility with PETSc
                topology->cells -= 1;
            } // if
            parser.ignore('}');
        } else if (0 == strcasecmp(token.c_str(), "material-ids")) {
            if (0 == topology->numCells) {
                const char* msg =
                    "Token 'count' must precede 'material-ids'.";
                throw std::runtime_error(msg);
            } // if
            materialIds->resize(topology->numCells);
            int label = 0;
            for (int iCell = 0; iCell < topology->numCells; ++iCell) {
                buffer.str(parser.next());
                buffer.clear();
                buffer >> label;
                buffer >> (*materialIds)[iCell];
            } // for

            // Zero index does NOT apply to materialIds.

            parser.ignore('}');
        } else {
            std::ostringstream msg;
            msg << "Could not parse '" << token << "' into an cells setting.";
            throw std::runtime_error(msg.str());
        } // else
        buffer.str(parser.next());
        buffer.clear();
        buffer >> token;
    } // while
    if (token != "}") {
        throw std::runtime_error("I/O error while parsing cells.");
    }
    topology->cellShape = pylith::meshio::MeshBuilder::cellShapeFromCorners(topology->dimension, topology->numCorners);

    // If no materials given, assign each cell material identifier of 0
    if ((0 == materialIds->size()) && (topology->numCells > 0)) {
        const int size = topology->numCells;
        materialIds->resize(size);
        (*materialIds) = 0;
    } // if

    PYLITH_METHOD_END;
} // readCells


// ------------------------------------------------------------------------------------------------
// Write mesh cells.
void
pylith::meshio::_MeshIOAscii::writeCells(std::ostream& fileout,
                                         const pylith::meshio::MeshBuilder::Topology& topology,
                                         const int_array& materialIds) {
    PYLITH_METHOD_BEGIN;

    fileout
        << "  cells = {\n"
        << "    count = " << topology.numCells << "\n"
        << "    num-corners = " << topology.numCorners << "\n"
        << "    simplices = {\n";

    for (int iCell = 0, i = 0; iCell < topology.numCells; ++iCell) {
        fileout << "      " << std::setw(8) << iCell;
        for (int iCorner = 0; iCorner < topology.numCorners; ++iCorner) {
            fileout << std::setw(8) << topology.cells[i++];
        }
        fileout << "\n";
    } // for
    fileout << "    }\n";

    // Write material identifiers
    assert(size_t(topology.numCells) == materialIds.size());
    fileout << "    material-ids = {\n";
    for (int iCell = 0; iCell < topology.numCells; ++iCell) {
        fileout << "      " << std::setw(8) << iCell;
        fileout << std::setw(4) << materialIds[iCell] << "\n";
    } // for
    fileout << "    }\n";

    fileout << "  }\n";

    PYLITH_METHOD_END;
} // writeCells


// ------------------------------------------------------------------------------------------------
// Read group of vertices.
void
pylith::meshio::_MeshIOAscii::readVertexGroup(int_array* points,
                                              std::string* name,
                                              spatialdata::utils::LineParser& parser,
                                              const bool useIndexZero) {
    PYLITH_METHOD_BEGIN;

    assert(points);
    assert(name);

    std::string token;
    std::istringstream buffer;
    const int maxIgnore = 1024;
    size_t groupSize = -1;
    buffer.str(parser.next());
    buffer.clear();
    buffer >> token;
    while (buffer.good() && token != "}") {
        if (0 == strcasecmp(token.c_str(), "name")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> std::ws;
            char cbuffer[maxIgnore];
            buffer.get(cbuffer, maxIgnore, '\n');
            *name = cbuffer;
        } else if (0 == strcasecmp(token.c_str(), "count")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> groupSize;
        } else if (0 == strcasecmp(token.c_str(), "indices")) {
            if (-1 == groupSize) {
                std::ostringstream msg;
                msg << "Tokens 'count' must precede 'indices'.";
                throw std::runtime_error(msg.str());
            } // if
            points->resize(groupSize);
            buffer.str(parser.next());
            buffer.clear();
            int i = 0;
            while (buffer.good() && i < groupSize) {
                buffer >> (*points)[i++];
                buffer >> std::ws;
                if (!buffer.good() && (i < groupSize)) {
                    buffer.str(parser.next());
                    buffer.clear();
                } // if
            } // while
            parser.ignore('}');
        } else {
            std::ostringstream msg;
            msg << "Could not parse '" << token << "' into a group setting.";
            throw std::runtime_error(msg.str());
        } // else
        buffer.str(parser.next());
        buffer.clear();
        buffer >> token;
    } // while
    if (token != "}") {
        std::ostringstream msg;
        msg << "I/O error while parsing group '" << *name << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!useIndexZero) {
        *points -= 1;
    } // if

    PYLITH_METHOD_END;
} // readVertexGroup


// ------------------------------------------------------------------------------------------------
// Write group of faces.
void
pylith::meshio::_MeshIOAscii::writeVertexGroup(std::ostream& fileout,
                                               const int_array& points,
                                               const char* name) {
    PYLITH_METHOD_BEGIN;

    fileout
        << "  vertex-group = {\n"
        << "    name = " << name << "\n"
        << "    count = " << points.size() << "\n"
        << "    indices = {\n";
    for (int iPoint = 0; iPoint < points.size(); ++iPoint) {
        fileout << "    " << points[iPoint] << "\n";
    } // for

    fileout
        << "    }\n"
        << "  }\n";

    PYLITH_METHOD_END;
} // writeVertexGroup


// ------------------------------------------------------------------------------------------------
// Read group of faces.
void
pylith::meshio::_MeshIOAscii::readFaceGroup(int_array* faceValues,
                                            std::string* name,
                                            spatialdata::utils::LineParser& parser,
                                            const pylith::meshio::MeshBuilder::shape_t faceShape,
                                            const bool useIndexZero) {
    PYLITH_METHOD_BEGIN;

    assert(faceValues);
    assert(name);

    const size_t numFaceValues = 1 + pylith::meshio::MeshBuilder::getNumVerticesFace(faceShape);

    std::string token;
    std::istringstream buffer;
    const int maxIgnore = 1024;
    size_t numFaces = -1;
    buffer.str(parser.next());
    buffer.clear();
    buffer >> token;
    while (buffer.good() && token != "}") {
        if (0 == strcasecmp(token.c_str(), "name")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> std::ws;
            char cbuffer[maxIgnore];
            buffer.get(cbuffer, maxIgnore, '\n');
            *name = cbuffer;
        } else if (0 == strcasecmp(token.c_str(), "count")) {
            buffer.ignore(maxIgnore, '=');
            buffer >> numFaces;
        } else if (0 == strcasecmp(token.c_str(), "indices")) {
            if (-1 == numFaces) {
                std::ostringstream msg;
                msg << "Tokens 'count' must precede 'indices'.";
                throw std::runtime_error(msg.str());
            } // if
            const size_t size = numFaces * numFaceValues;
            faceValues->resize(size);
            buffer.str(parser.next());
            buffer.clear();
            int i = 0;
            while (buffer.good() && i < size) {
                buffer >> (*faceValues)[i++];
                buffer >> std::ws;
                if (!buffer.good() && (i < size)) {
                    buffer.str(parser.next());
                    buffer.clear();
                } // if
            } // while
            parser.ignore('}');
        } else {
            std::ostringstream msg;
            msg << "Could not parse '" << token << "' into a group setting.";
            throw std::runtime_error(msg.str());
        } // else
        buffer.str(parser.next());
        buffer.clear();
        buffer >> token;
    } // while
    if (token != "}") {
        std::ostringstream msg;
        msg << "I/O error while parsing group '" << *name << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!useIndexZero) {
        *faceValues -= 1;
    } // if

    PYLITH_METHOD_END;
} // readFaceGroup


// ------------------------------------------------------------------------------------------------
// Write group of faces.
void
pylith::meshio::_MeshIOAscii::writeFaceGroup(std::ostream& fileout,
                                             const int_array& faceValues,
                                             const pylith::meshio::MeshBuilder::shape_t faceShape,
                                             const char* name) {
    PYLITH_METHOD_BEGIN;

    const size_t numFaceValues = 1 + pylith::meshio::MeshBuilder::getNumVerticesFace(faceShape);
    assert(0 == faceValues.size() % numFaceValues);
    const size_t numFaces = faceValues.size() / numFaceValues;

    fileout
        << "  face-group = {\n"
        << "    name = " << name << "\n"
        << "    count = " << numFaces << "\n"
        << "    indices = {\n";
    for (int iFace = 0; iFace < numFaces; ++iFace) {
        fileout << "    ";
        for (int iValue = 0; iValue < numFaceValues; ++iValue) {
            fileout << "  " << faceValues[numFaceValues*iFace+iValue];
        } // for
        fileout << "\n";
    } // for

    fileout
        << "    }\n"
        << "  }\n";

    PYLITH_METHOD_END;
} // writeFaceGroup


// End of file
