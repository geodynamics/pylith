// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MeshIOCubit.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder
#include "ExodusII.hh" // USES ExodusII

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "petsc.h" // USES MPI_Comm

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOCubit::MeshIOCubit(void) :
    _filename(""),
    _useNodesetNames(true) { // constructor
    PyreComponent::setName("meshiocubit");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOCubit::~MeshIOCubit(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOCubit::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOCubit::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");

    assert(_mesh);

    const int commRank = _mesh->getCommRank();
    int meshDim = 0;
    int spaceDim = 0;
    int numVertices = 0;
    int numCells = 0;
    int numCorners = 0;
    scalar_array coordinates;
    int_array cells;
    int_array materialIds;
    PetscErrorCode err = 0;

    if (0 == commRank) {
        try {
            ExodusII exofile(_filename.c_str());

            const int meshDim = exofile.getDim("num_dim");

            _readVertices(exofile, &coordinates, &numVertices, &spaceDim);
            err = MPI_Bcast(&spaceDim, 1, MPI_INT, 0, _mesh->getComm());PYLITH_CHECK_ERROR(err);
            _readCells(exofile, &cells, &materialIds, &numCells, &numCorners);
            _orientCells(&cells, numCells, numCorners, meshDim);
            MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners, meshDim);
            _setMaterials(materialIds);

            _readGroups(exofile);
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
        err = MPI_Bcast(&spaceDim, 1, MPI_INT, 0, _mesh->getComm());PYLITH_CHECK_ERROR(err);
        MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners, meshDim);
        _setMaterials(materialIds);
    }
    _distributeGroups();

    PYLITH_METHOD_END;
} // read


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOCubit::_write(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_write()");

    ExodusII exofile(_filename.c_str());

    _writeDimensions(exofile);
    _writeVariables(exofile);
    _writeAttributes(exofile);

    PYLITH_METHOD_END;
} // write


// ---------------------------------------------------------------------------------------------------------------------
// Read mesh vertices.
void
pylith::meshio::MeshIOCubit::_readVertices(ExodusII& exofile,
                                           scalar_array* coordinates,
                                           int* numVertices,
                                           int* numDims) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_readVertices(exofile="<<typeid(exofile).name()<<", coordinates="<<coordinates<<", numVertices="<<numVertices<<", numDims="<<numDims<<")");

    assert(coordinates);
    assert(numVertices);
    assert(numDims);

    // Space dimension
    *numDims = exofile.getDim("num_dim");

    // Number of vertices
    *numVertices = exofile.getDim("num_nodes");

    PYLITH_COMPONENT_INFO("Reading " << *numVertices << " vertices.");

    if (exofile.hasVar("coord", NULL)) {
        const int ndims = 2;
        int dims[2];
        dims[0] = *numDims;
        dims[1] = *numVertices;
        scalar_array buffer(*numVertices * *numDims);
        exofile.getVar(&buffer[0], dims, ndims, "coord");

        coordinates->resize(*numVertices * *numDims);
        for (int iVertex = 0; iVertex < *numVertices; ++iVertex) {
            for (int iDim = 0; iDim < *numDims; ++iDim) {
                (*coordinates)[iVertex*(*numDims)+iDim] =
                    buffer[iDim*(*numVertices)+iVertex];
            }
        }

    } else {
        const char* coordNames[3] = { "coordx", "coordy", "coordz" };

        coordinates->resize(*numVertices * *numDims);
        scalar_array buffer(*numVertices);

        const int ndims = 1;
        int dims[1];
        dims[0] = *numVertices;

        for (int i = 0; i < *numDims; ++i) {
            exofile.getVar(&buffer[0], dims, ndims, coordNames[i]);

            for (int iVertex = 0; iVertex < *numVertices; ++iVertex) {
                (*coordinates)[iVertex*(*numDims)+i] = buffer[iVertex];
            }
        } // for
    } // else

    PYLITH_METHOD_END;
} // _readVertices


// ---------------------------------------------------------------------------------------------------------------------
// Read mesh cells.
void
pylith::meshio::MeshIOCubit::_readCells(ExodusII& exofile,
                                        int_array* cells,
                                        int_array* materialIds,
                                        int* numCells,
                                        int* numCorners) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_readCells(exofile="<<typeid(exofile).name()<<", cells="<<cells<<", materialIds="<<materialIds<<", numCells="<<numCells<<", numCorners="<<numCorners<<")");

    assert(cells);
    assert(materialIds);
    assert(numCells);
    assert(numCorners);

    *numCells = exofile.getDim("num_elem");
    const int numMaterials = exofile.getDim("num_el_blk");

    PYLITH_COMPONENT_INFO("Reading " << *numCells << " cells in " << numMaterials << " blocks.");

    int_array blockIds(numMaterials);
    int ndims = 1;
    int dims[2];
    dims[0] = numMaterials;
    dims[1] = 0;
    exofile.getVar(&blockIds[0], dims, ndims, "eb_prop1");

    materialIds->resize(*numCells);
    *numCorners = 0;
    for (int iMaterial = 0, index = 0; iMaterial < numMaterials; ++iMaterial) {
        std::ostringstream varname;
        varname << "num_nod_per_el" << iMaterial+1;
        if (0 == *numCorners) {
            *numCorners = exofile.getDim(varname.str().c_str());
            const int size = (*numCells) * (*numCorners);
            cells->resize(size);
        } else if (exofile.getDim(varname.str().c_str()) != *numCorners) {
            std::ostringstream msg;
            msg << "All materials must have the same number of vertices per cell.\n"
                << "Expected " << *numCorners << " vertices per cell, but block "
                << blockIds[iMaterial] << " has "
                << exofile.getDim(varname.str().c_str())
                << " vertices.";
            throw std::runtime_error(msg.str());
        } // if

        varname.str("");
        varname << "num_el_in_blk" << iMaterial+1;
        const int blockSize = exofile.getDim(varname.str().c_str());

        varname.str("");
        varname << "connect" << iMaterial+1;
        ndims = 2;
        dims[0] = blockSize;
        dims[1] = *numCorners;
        exofile.getVar(&(*cells)[index* (*numCorners)], dims, ndims,
                       varname.str().c_str());

        for (int i = 0; i < blockSize; ++i) {
            (*materialIds)[index+i] = blockIds[iMaterial];
        }

        index += blockSize;
    } // for

    *cells -= 1; // use zero index

    PYLITH_METHOD_END;
} // _readCells


// ---------------------------------------------------------------------------------------------------------------------
// Read mesh groups.
void
pylith::meshio::MeshIOCubit::_readGroups(ExodusII& exofile) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_readGroups(exofile="<<typeid(exofile).name()<<")");

    const int numGroups = exofile.getDim("num_node_sets");

    PYLITH_COMPONENT_INFO("Found " << numGroups << " node sets.");

    int_array ids(numGroups);
    int ndims = 1;
    int dims[2];
    dims[0] = numGroups;
    dims[1] = 0;
    exofile.getVar(&ids[0], dims, ndims, "ns_prop1");

    string_vector groupNames(numGroups);

    if (_useNodesetNames) {
        exofile.getVar(&groupNames, numGroups, "ns_names");
    } // if

    for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
        std::ostringstream varname;
        varname << "num_nod_ns" << iGroup+1;
        const size_t nodesetSize = exofile.getDim(varname.str().c_str());
        int_array points(nodesetSize);

        varname.str("");
        varname << "node_ns" << iGroup+1;
        ndims = 1;
        dims[0] = nodesetSize;

        PYLITH_COMPONENT_INFO("Reading node set '" << groupNames[iGroup] << "' with id " << ids[iGroup] << " containing " << nodesetSize << " nodes.");
        exofile.getVar(&points[0], dims, ndims, varname.str().c_str());

        std::sort(&points[0], &points[0]+nodesetSize);
        points -= 1; // use zero index

        GroupPtType type = VERTEX;
        if (_useNodesetNames) {
            _setGroup(groupNames[iGroup], type, points);
        } else {
            std::ostringstream name;
            name << ids[iGroup];
            _setGroup(name.str().c_str(), type, points);
        } // if/else
    } // for

    PYLITH_METHOD_END;
} // _readGroups


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh dimensions.
void
pylith::meshio::MeshIOCubit::_writeDimensions(ExodusII& exofile) const {
    throw std::logic_error("MeshIOCubit::_writeDimensions() not implemented.");
} // _writeDimensions


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh variables.
void
pylith::meshio::MeshIOCubit::_writeVariables(ExodusII& exofile) const {
    throw std::logic_error("MeshIOCubit::_writeVariables() not implemented.");
} // _writeVariables


// ---------------------------------------------------------------------------------------------------------------------
// Write mesh attributes.
void
pylith::meshio::MeshIOCubit::_writeAttributes(ExodusII& exofile) const {
    throw std::logic_error("MeshIOCubit::_writeAttributes() not implemented.");
} // _writeAttributes


// ---------------------------------------------------------------------------------------------------------------------
// Reorder vertices in cells to match PyLith conventions.
void
pylith::meshio::MeshIOCubit::_orientCells(int_array* const cells,
                                          const int numCells,
                                          const int numCorners,
                                          const int meshDim) {
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(cells->size() == size_t(numCells*numCorners));

    if ((2 == meshDim) && (4 == numCorners)) { // QUAD4
        // do nothing

    } else if ((3 == meshDim) && (8 == numCorners)) { // HEX8
        // do nothing

    } else if ((2 == meshDim) && (6 == numCorners)) { // TRI6
        // CUBIT
        // corners,
        // bottom edges, middle edges, top edges

        // PyLith
        // bottom edge, right edge, left edge, corners

        // Permutation: 3, 4, 5, 0, 1, 2
        int tmp = 0;
        for (int iCell = 0; iCell < numCells; ++iCell) {
            const int ii = iCell*numCorners;
            tmp = (*cells)[ii+0];
            (*cells)[ii+0] = (*cells)[ii+3];
            (*cells)[ii+3] = tmp;

            tmp = (*cells)[ii+1];
            (*cells)[ii+1] = (*cells)[ii+4];
            (*cells)[ii+4] = tmp;

            tmp = (*cells)[ii+2];
            (*cells)[ii+2] = (*cells)[ii+5];
            (*cells)[ii+5] = tmp;
        } // for

    } else if ((3 == meshDim) && (27 == numCorners)) { // HEX27
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
        for (int iCell = 0; iCell < numCells; ++iCell) {
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

            tmp = (*cells)[i12];
            (*cells)[i12] = (*cells)[i16];
            (*cells)[i16] = tmp;

            tmp = (*cells)[i13];
            (*cells)[i13] = (*cells)[i17];
            (*cells)[i17] = tmp;

            tmp = (*cells)[i14];
            (*cells)[i14] = (*cells)[i18];
            (*cells)[i18] = tmp;

            tmp = (*cells)[i15];
            (*cells)[i15] = (*cells)[i19];
            (*cells)[i19] = tmp;

            tmp = (*cells)[i20];
            (*cells)[i20] = (*cells)[i23];
            (*cells)[i23] = (*cells)[i26];
            (*cells)[i26] = tmp;

            tmp = (*cells)[i21];
            (*cells)[i21] = (*cells)[i24];
            (*cells)[i24] = tmp;

            tmp = (*cells)[i22];
            (*cells)[i22] = (*cells)[i25];
            (*cells)[i25] = tmp;
        } // for
    } // if/else

    PYLITH_METHOD_END;
} // _orientCells


// End of file
