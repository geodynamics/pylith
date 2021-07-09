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

#include "MeshIOLagrit.hh" // implementation of class methods

#include "GMVFileAscii.hh" // USES GMVFileAscii
#include "GMVFileBinary.hh" // USES GMVFileBinary
#include "PsetFileAscii.hh" // USES PsetFileAscii
#include "PsetFileBinary.hh" // USES PsetFileBinary
#include "MeshBuilder.hh" // USES MeshBuilder

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/array.hh" // USES scalar_array, int_array

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error()

// ---------------------------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOLagrit::MeshIOLagrit(void) :
    _filenameGmv(""),
    _filenamePset(""),
    _flipEndian(false),
    _ioInt32(true),
    _isRecordHeader32Bit(true) {
    PyreComponent::setName("meshiolagrit");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOLagrit::~MeshIOLagrit(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOLagrit::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Unpickle mesh
void
pylith::meshio::MeshIOLagrit::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");

    const int commRank = _mesh->getCommRank();
    int meshDim = 0;
    int spaceDim = 0;
    int numVertices = 0;
    int numCells = 0;
    int numCorners = 0;
    scalar_array coordinates;
    int_array cells;
    int_array materialIds;

    if (!commRank) {
        if (GMVFile::isAscii(_filenameGmv.c_str())) {
            GMVFileAscii filein(_filenameGmv.c_str());
            filein.read(&coordinates, &cells, &materialIds,
                        &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
            _orientCellsAscii(&cells, numCells, numCorners, meshDim);
        } else {
            GMVFileBinary filein(_filenameGmv.c_str(), _flipEndian);
            filein.read(&coordinates, &cells, &materialIds,
                        &meshDim, &spaceDim, &numVertices, &numCells, &numCorners);
            _orientCellsBinary(&cells, numCells, numCorners, meshDim);
        } // if/else
    }
    MeshBuilder::buildMesh(_mesh, &coordinates, numVertices, spaceDim, cells, numCells, numCorners, meshDim);
    _setMaterials(materialIds);

    if (0 == commRank) {
        std::vector<PsetFile::Pset> groups;
        if (PsetFile::isAscii(_filenamePset.c_str())) {
            PsetFileAscii filein(_filenamePset.c_str());
            filein.read(&groups);
        } else {
            PsetFileBinary filein(_filenamePset.c_str(),
                                  _flipEndian,
                                  _ioInt32,
                                  _isRecordHeader32Bit);
            filein.read(&groups);
        } // if/else
        GroupPtType type = VERTEX;
        const int numGroups = groups.size();
        for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
            _setGroup(groups[iGroup].name, type, groups[iGroup].points);
        }
    }
    _distributeGroups();

    PYLITH_METHOD_END;
} // _read


// ---------------------------------------------------------------------------------------------------------------------
// Pickle mesh
void
pylith::meshio::MeshIOLagrit::_write(void) const {
    throw std::logic_error("MeshIOLagrit::_write not implemented.");
} // _write


// ---------------------------------------------------------------------------------------------------------------------
// Reorder vertices in cells from ASCII GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsAscii(int_array* const cells,
                                                const int numCells,
                                                const int numCorners,
                                                const int meshDim) {
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(cells->size() == size_t(numCells*numCorners));

    if ((3 == meshDim) && (4 == numCorners)) { // TET
        for (int iCell = 0; iCell < numCells; ++iCell) {
            const int i1 = iCell*numCorners+1;
            const int i2 = iCell*numCorners+2;
            const int tmp = (*cells)[i1];
            (*cells)[i1] = (*cells)[i2];
            (*cells)[i2] = tmp;
        } // for

    } // if
    PYLITH_METHOD_END;
} // _orientCellsAscii


// ---------------------------------------------------------------------------------------------------------------------
// Reorder vertices in cells from binary GMV file to match PyLith
// conventions.
void
pylith::meshio::MeshIOLagrit::_orientCellsBinary(int_array* const cells,
                                                 const int numCells,
                                                 const int numCorners,
                                                 const int meshDim) {
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(cells->size() == size_t(numCells*numCorners));

    if ((3 == meshDim) && (4 == numCorners)) { // TET
        // do nothing
    } // if
    PYLITH_METHOD_END;
} // _orientCellsBinary


// End of file
