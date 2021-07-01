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

#include "PsetFileAscii.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <cstring> // USES strcmp()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception

// ----------------------------------------------------------------------
const char* pylith::meshio::PsetFileAscii::_HEADER = "pset ascii";

// ----------------------------------------------------------------------
// Constructor with name of Pset file.
pylith::meshio::PsetFileAscii::PsetFileAscii(const char* filename) :
    PsetFile(filename) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Default destructor
pylith::meshio::PsetFileAscii::~PsetFileAscii(void) { // destructor
} // destructor


// ----------------------------------------------------------------------
// Read ASCII Pset file.
void
pylith::meshio::PsetFileAscii::read(std::vector<Pset>* groups) { // read
    PYLITH_METHOD_BEGIN;

    assert(groups);

    pythia::journal::info_t info("meshio");

    std::ifstream fin(_filename.c_str(), std::ios::in);
    if (!(fin.is_open() && fin.good())) {
        std::ostringstream msg;
        msg
            << "Could not open ASCII Pset file '" << _filename
            << "' for reading.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Reading ASCII Pset file '" << _filename << "'." << pythia::journal::endl;

    _readHeader(fin);

    // Read number of psets
    int numGroups = 0;
    fin >> numGroups;
    groups->resize(numGroups);

    // Read groups
    info << pythia::journal::at(__HERE__)
         << "Reading " << numGroups << " point sets from file." << pythia::journal::endl;
    for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
        _readPset(fin, &(*groups)[iGroup]);
    }

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write ASCII Pset file.
void
pylith::meshio::PsetFileAscii::write(const std::vector<Pset>& groups) { // write
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");

    std::ofstream fout(_filename.c_str(), std::ios::out);
    if (!(fout.is_open() && fout.good())) {
        std::ostringstream msg;
        msg
            << "Could not open ASCII Pset file '" << _filename
            << "' for writing.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Writing ASCII Pset file '" << _filename << "'." << pythia::journal::endl;

    _writeHeader(fout);

    // Write number of groups
    const int numGroups = groups.size();
    fout << std::setw(4) << numGroups << std::endl;

    // Write groups
    info << pythia::journal::at(__HERE__)
         << "Writing " << numGroups << " point sets to file." << pythia::journal::endl;
    for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
        _writePset(fout, groups[iGroup]);
    }

    PYLITH_METHOD_END;
} // write


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_readHeader(std::ifstream& fin) { // _readHeader
    PYLITH_METHOD_BEGIN;

    const int headerLen = strlen(_HEADER)+1;
    char buffer[headerLen];
    fin.get(buffer, headerLen, '\n');
    if (strcmp(_HEADER, buffer)) {
        std::ostringstream msg;
        msg
            << "Header in ASCII Pset file '" << buffer
            << "' does not match anticipated header '"
            << _HEADER << "'";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // _readHeader


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_writeHeader(std::ofstream& fout) { // _writeHeader
    PYLITH_METHOD_BEGIN;

    fout << _HEADER << std::endl;

    PYLITH_METHOD_END;
} // _writeHeader


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_readPset(std::ifstream& fin,
                                         Pset* group) { // _readPset
    PYLITH_METHOD_BEGIN;

    assert(group);

    pythia::journal::info_t info("meshio");

    int size = 0;
    fin >> group->name >> group->id >> size;
    info << pythia::journal::at(__HERE__)
         << "Reading point set '" << group->name << "' with " << size
         << " points." << pythia::journal::endl;

    group->points.resize(size);
    for (int i = 0; i < size; ++i) {
        fin >> group->points[i];
    }

    group->points -= 1; // use zero base

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readPset


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileAscii::_writePset(std::ofstream& fout,
                                          const Pset& group) { // _writePset
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");
    const int size = group.points.size();

    info << pythia::journal::at(__HERE__)
         << "Writing point set '" << group.name << "' with " << size
         << " points." << pythia::journal::endl;

    fout << group.name << " " << group.id << "  " << size << std::endl;

    const int numCols = 10;
    for (int i = 0, iCol = 0; i < size; ++i) {
        fout << std::setw(8) << 1+group.points[i];
        if (++iCol == numCols) {
            fout << std::endl;
            iCol = 0;
        } // if
    } // for

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _writePset


// End of file
