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

#include "GMVFileAscii.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES pythia::journal

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <cstring> // USES strcmp()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::GMVFileAscii::_HEADER = "gmvinput ascii";

// ----------------------------------------------------------------------
// Constructor with name of GMV file.
pylith::meshio::GMVFileAscii::GMVFileAscii(const char* filename) :
    GMVFile(filename) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Default destructor
pylith::meshio::GMVFileAscii::~GMVFileAscii(void) { // destructor
} // destructor


// ----------------------------------------------------------------------
// Read ASCII GMV file.
void
pylith::meshio::GMVFileAscii::read(scalar_array* coordinates,
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

    pythia::journal::info_t info("gmvfile");

    std::ifstream fin(_filename.c_str(), std::ios::in);
    if (!(fin.is_open() && fin.good())) {
        std::ostringstream msg;
        msg
            << "Could not open ASCII GMV file '" << _filename
            << "' for reading.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Reading ASCII GMV file '" << _filename << "'." << pythia::journal::endl;

    _readHeader(fin);

    std::string token;
    while (fin >> token) {
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
    } // while

    assert(coordinates->size() == size_t((*numVertices) * (*spaceDim)));
    assert(cells->size() == size_t((*numCells) * (*numCorners)));
    assert(materialIds->size() == size_t(*numCells));

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write ASCII GMV file.
void
pylith::meshio::GMVFileAscii::write(const scalar_array& coordinates,
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
pylith::meshio::GMVFileAscii::_readHeader(std::ifstream& fin) { // _readHeader
    PYLITH_METHOD_BEGIN;

    const int headerLen = strlen(_HEADER)+1;
    char buffer[headerLen];
    fin.get(buffer, headerLen, '\n');
    if (strcmp(_HEADER, buffer)) {
        std::ostringstream msg;
        msg
            << "Header in ASCII GMV file '" << buffer
            << "' does not match anticipated header '"
            << _HEADER << "'";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // _readHeader


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readVertices(std::ifstream& fin,
                                            scalar_array* coordinates,
                                            int* numVertices,
                                            int* spaceDim) { // _readVertices
    PYLITH_METHOD_BEGIN;

    assert(coordinates);
    assert(numVertices);
    assert(spaceDim);

    *spaceDim = 3;

    pythia::journal::info_t info("gmvfile");

    fin >> *numVertices;
    info << pythia::journal::at(__HERE__)
         << "Reading " << *numVertices << " nodes." << pythia::journal::endl;

    coordinates->resize(*numVertices * (*spaceDim));
    // NOTE: Order of loops is different than what we usually have
    for (int iDim = 0; iDim < *spaceDim; ++iDim) {
        for (int iVertex = 0; iVertex < *numVertices; ++iVertex) {
            fin >> (*coordinates)[iVertex*(*spaceDim)+iDim];
        }
    }

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // readVertices


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readCells(std::ifstream& fin,
                                         int_array* cells,
                                         int* numCells,
                                         int* numCorners) { // readCells
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(numCells);
    assert(numCorners);

    pythia::journal::info_t info("gmvfile");

    *numCorners = 0;

    fin >> *numCells;
    std::string cellString = "";
    info << pythia::journal::at(__HERE__)
         << "Reading " << numCells << " cells." << pythia::journal::endl;
    for (int iCell = 0; iCell < *numCells; ++iCell) {
        std::string cellStringCur;
        int numCornersCur = 0;
        fin >> cellStringCur;
        fin >> numCornersCur;
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
        for (int iCorner = 0; iCorner < *numCorners; ++iCorner) {
            fin >> (*cells)[iCell*(*numCorners)+iCorner];
        }
        fin >> std::ws;
    } // for

    *cells -= 1; // use zero base

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // readCells


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readVariables(std::ifstream& fin,
                                             const int numVertices,
                                             const int numCells) { // _readVariables
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("gmvfile");

    info << pythia::journal::at(__HERE__)
         << "Reading variables..." << pythia::journal::endl;

    std::string varName;
    fin >> varName;
    while ("endvars" != varName && !fin.eof() && fin.good()) {
        int varType = 0;
        fin >> varType;
        if (1 == varType) { // variables/attributes associated with vertices
            const int numVars = 1;
            scalar_array vals(numVertices*numVars);
            for (int iVertex = 0; iVertex < numVertices; ++iVertex) {
                fin >> vals[iVertex];
            }
        } else { // variables/attributes associated with cells
            const int numVars = 1;
            scalar_array vals(numCells*numVars);
            for (int iCell = 0; iCell < numCells; ++iCell) {
                fin >> vals[iCell];
            }
        } // else
        fin >> varName;
    } // while

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readVariables


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readFlags(std::ifstream& fin,
                                         const int numVertices,
                                         const int numCells) { // _readFlags
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("gmvfile");

    info << pythia::journal::at(__HERE__)
         << "Reading flags..." << pythia::journal::endl;

    std::string varName;
    fin >> varName;
    while ("endflag" != varName && !fin.eof() && fin.good()) {
        int varType = 0;
        int numFlags = 0;
        fin >> numFlags >> varType;
        for (int iFlag = 0; iFlag < numFlags; ++iFlag) {
            std::string flagName;
            fin >> flagName;
        } // for
        if (1 == varType) { // flag associated with vertices
            int flag;
            for (int iVertex = 0; iVertex < numVertices; ++iVertex) {
                fin >> flag;
            }
        } else { // flag associated with cells
            int flag;
            for (int iCell = 0; iCell < numCells; ++iCell) {
                fin >> flag;
            }
        } // else
        fin >> varName;
    } // while

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readFlags


// ----------------------------------------------------------------------
void
pylith::meshio::GMVFileAscii::_readMaterials(std::ifstream& fin,
                                             int_array* materialIds,
                                             const int numVertices,
                                             const int numCells) { // _readMaterials
    PYLITH_METHOD_BEGIN;

    assert(materialIds);

    pythia::journal::info_t info("gmvfile");
    info << pythia::journal::at(__HERE__)
         << "Reading materials..." << pythia::journal::endl;

    int numMaterials = 0;
    int dataType = 0;
    fin >> numMaterials >> dataType;

    std::string name;
    for (int iMat = 0; iMat < numMaterials; ++iMat) {
        fin >> name;
    }

    if (0 == dataType) { // material associated with elements
        materialIds->resize(numCells);
        for (int iCell = 0; iCell < numCells; ++iCell) {
            fin >> (*materialIds)[iCell];
        }
    } else { // material associated with nodes
        int_array materials(numVertices);
        for (int iVertex = 0; iVertex < numVertices; ++iVertex) {
            fin >> materials[iVertex];
        }
    } // else

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readMaterials


// End of file
