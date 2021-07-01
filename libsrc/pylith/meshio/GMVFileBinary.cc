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

#include "GMVFileBinary.hh" // implementation of class methods

#include "BinaryIO.hh" // USES readString()

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh"

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::GMVFileBinary::_HEADER = "gmvinputieee    ";

// ----------------------------------------------------------------------
// Constructor with name of GMV file.
pylith::meshio::GMVFileBinary::GMVFileBinary(const char* filename,
                                             const bool flipEndian) :
    GMVFile(filename),
    _flipEndian(flipEndian) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Default destructor
pylith::meshio::GMVFileBinary::~GMVFileBinary(void) { // destructor
} // destructor


// ----------------------------------------------------------------------
// Read binary GMV file.
void
pylith::meshio::GMVFileBinary::read(scalar_array* coordinates,
                                    int_array* cells,
                                    int_array* materialIds,
                                    int* meshDim,
                                    int* spaceDim,
                                    int* numVertices,
                                    int* numCells,
                                    int* numCorners) { // read
    PYLITH_METHOD_BEGIN;

    assert(coordinates);
    assert(cells);
    assert(materialIds);
    assert(meshDim);
    assert(spaceDim);
    assert(numVertices);
    assert(numCells);
    assert(numCorners);

    *meshDim = 3;

    pythia::journal::info_t info("meshio");

    std::ifstream fin(_filename.c_str(), std::ios::in | std::ios::binary);
    if (!(fin.is_open() && fin.good())) {
        std::ostringstream msg;
        msg
            << "Could not open binary GMV file '" << _filename
            << "' for reading.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Reading binary GMV file '" << _filename << "'." << pythia::journal::endl;

    _readHeader(fin);

    const int tokenLen = 8;
    std::string token = BinaryIO::readString(fin, tokenLen);
    while (!fin.eof() && fin.good()) {
        if (token == "nodes") {
            _readVertices(fin, coordinates, numVertices, spaceDim);
        } else if (token == "cells") {
            _readCells(fin, cells, numCells, numCorners);
        } else if (token == "variable") {
            _readVariables(fin, *numVertices, *numCells);
        } else if (token == "flags") {
            _readFlags(fin, *numVertices, *numCells);
        } else if (token == "material") {
            _readMaterials(fin, materialIds, *numVertices, *numCells);
        }
        token = BinaryIO::readString(fin, tokenLen);
    } // while

    assert(coordinates->size() == size_t((*numVertices) * (*spaceDim)));
    assert(cells->size() == size_t((*numCells) * (*numCorners)));
    assert(materialIds->size() == size_t(*numCells));

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write binary GMV file.
void
pylith::meshio::GMVFileBinary::write(const scalar_array& coordinates,
                                     const int_array& cells,
                                     const int_array& materialIds,
                                     const int meshDim,
                                     const int spaceDim,
                                     const int numVertices,
                                     const int numCells,
                                     const int numCorners) { // write
    PYLITH_METHOD_BEGIN;

    assert(coordinates.size() == size_t(numVertices * spaceDim));
    assert(cells.size() == size_t(numCells * numCorners));
    assert(materialIds.size() == size_t(numCells));

#if 0 // NOT YET IMPLEMENTED
    _writeHeader();
    _writeVertices(coordinates);
    _writeCells(cells);
    _writeMaterials(materialIds);
#endif

    PYLITH_METHOD_END;
} // write


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readHeader(std::ifstream& fin) { // _readHeader
    PYLITH_METHOD_BEGIN;

    std::string header = BinaryIO::readString(fin, strlen(_HEADER));
    std::string headerE = _HEADER;
    headerE = headerE.substr(0, headerE.find_first_of(" "));
    if (headerE != header) {
        std::ostringstream msg;
        msg
            << "Header in binary GMV file '" << header
            << "' does not match anticipated header '" << headerE << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // _readHeader


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readVertices(std::ifstream& fin,
                                             scalar_array* coordinates,
                                             int* numVertices,
                                             int* spaceDim) { // _readVertices
    PYLITH_METHOD_BEGIN;

    assert(coordinates);
    assert(numVertices);
    assert(spaceDim);

    *spaceDim = 3;

    pythia::journal::info_t info("meshio");

    fin.read((char*) numVertices, sizeof(int));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) numVertices, 1, sizeof(int));
    }
    assert(*numVertices > 0);
    info << pythia::journal::at(__HERE__)
         << "Reading " << *numVertices << " nodes." << pythia::journal::endl;

    const int size = (*numVertices) * (*spaceDim);
    float_array buffer(size);

    fin.read((char*) &buffer[0], size*sizeof(float));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &buffer[0], size, sizeof(float));
    }

    coordinates->resize(size);
    // switch from column major to row major order
    for (int iDim = 0, i = 0; iDim < *spaceDim; ++iDim) {
        for (int iVertex = 0; iVertex < *numVertices; ++iVertex) {
            (*coordinates)[iVertex*(*spaceDim)+iDim] = buffer[i++];
        }
    }

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readVertices


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readCells(std::ifstream& fin,
                                          int_array* cells,
                                          int* numCells,
                                          int* numCorners) { // _readCells
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(numCells);
    assert(numCorners);

    pythia::journal::info_t info("meshio");

    fin.read((char*) numCells, sizeof(int));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) numCells, 1, sizeof(int));
    }
    info << pythia::journal::at(__HERE__)
         << "Reading " << *numCells << " cells." << pythia::journal::endl;
    std::string cellString = "";
    for (int iCell = 0; iCell < *numCells; ++iCell) {
        const int stringLen = 8;
        char cellStringCur[stringLen+1];
        int numCornersCur = 0;
        fin.read((char*) cellStringCur, sizeof(char)*stringLen);
        cellStringCur[stringLen] = '\0';
        fin.read((char*) &numCornersCur, sizeof(int));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numCornersCur, 1, sizeof(int));
        }
        if (*numCorners) {
            if (cellStringCur != cellString) {
                std::ostringstream msg;
                msg
                    << "Mutiple element types not supported. Found element types '"
                    << cellString << "' and '" << cellStringCur << "' in GMV file '"
                    << _filename << "'.";
                throw std::runtime_error(msg.str());
            } // if
            assert(*numCorners == numCornersCur);
        } else {
            cellString = cellStringCur;
            *numCorners = numCornersCur;
            cells->resize((*numCells) * (*numCorners));
        } // if/else
        fin.read((char*) &(*cells)[iCell*numCornersCur],
                 sizeof(int)*numCornersCur);
    } // for
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &(*cells)[0],
                                (*numCells)*(*numCorners), sizeof(int));
    }

    *cells -= 1; // use zero base

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readCells


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readVariables(std::ifstream& fin,
                                              const int numVertices,
                                              const int numCells) { // _readVariables
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");
    info << pythia::journal::at(__HERE__)
         << "Reading variables..." << pythia::journal::endl;

    const int varNameLen = 8;
    std::string varName = BinaryIO::readString(fin, varNameLen);
    while ("endvars" != varName && !fin.eof() && fin.good()) {
        int varType = 0;
        fin.read((char*) &varType, sizeof(int));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &varType, 1, sizeof(varType));
        }
        if (1 == varType) { // variable/attribute associated with vertices
            float_array vals(numVertices);
            fin.read((char*) &vals[0], sizeof(float)*numVertices);
        } else { // variable/attribute associated with cells
            float_array vals(numCells);
            fin.read((char*) &vals[0], sizeof(float)*numCells);
        } // else
        varName = BinaryIO::readString(fin, varNameLen);
    } // while

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readVariables


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readFlags(std::ifstream& fin,
                                          const int numVertices,
                                          const int numCells) { // _readFlags
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");
    info << pythia::journal::at(__HERE__)
         << "Reading flags..." << pythia::journal::endl;

    const int varNameLen = 8;
    std::string varName = BinaryIO::readString(fin, varNameLen);
    while ("endflag" != varName && !fin.eof() && fin.good()) {
        int numFlags = 0;
        fin.read((char*) &numFlags, sizeof(int));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numFlags, 1, sizeof(numFlags));
        }
        int varType = 0;
        fin.read((char*) &varType, sizeof(int));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &varType, 1, sizeof(varType));
        }
        for (int iFlag = 0; iFlag < numFlags; ++iFlag) {
            const int flagNameLen = 8;
            std::string flagName = BinaryIO::readString(fin, flagNameLen);
        } // for
        if (1 == varType) { // flag associated with vertices
            int_array buffer(numVertices);
            fin.read((char*) &buffer[0], sizeof(int)*numVertices);
        } else { // flag associated with cells
            int_array buffer(numVertices);
            fin.read((char*) &buffer[0], sizeof(int)*numCells);
        } // else
        varName = BinaryIO::readString(fin, varNameLen);
    } // while

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readFlags


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileBinary::_readMaterials(std::ifstream& fin,
                                              int_array* materialIds,
                                              const int numVertices,
                                              const int numCells) { // _readMaterials
    PYLITH_METHOD_BEGIN;

    assert(materialIds);

    pythia::journal::info_t info("meshio");
    info << pythia::journal::at(__HERE__)
         << "Reading materials..." << pythia::journal::endl;

    int numMaterials = 0;
    fin.read((char*) &numMaterials, sizeof(int));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &numMaterials, 1, sizeof(numMaterials));
    }
    int dataType = 0;
    fin.read((char*) &dataType, sizeof(int));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &dataType, 1, sizeof(dataType));
    }

    for (int iMat = 0; iMat < numMaterials; ++iMat) {
        const int nameLen = 8;
        std::string name = BinaryIO::readString(fin, nameLen);
    } // for

    if (0 == dataType) { // materials associated with cells
        materialIds->resize(numCells);
        fin.read((char*) &(*materialIds)[0], sizeof(int)*numCells);
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &(*materialIds)[0], numCells,
                                    sizeof(int));
        }
    } else { // materials associated with vertices
        int_array buffer(numVertices);
        fin.read((char*) &buffer[0], sizeof(int)*numVertices);
    } // else

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readMaterials


// End of file
