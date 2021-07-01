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

#include "PsetFileBinary.hh" // implementation of class methods

#include "BinaryIO.hh" // USES readString()

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // pythia::journal

#include <fstream> // USES std::ifstream
#include <iomanip> // USES std::setw()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::exception
#include <sys/types.h> // USES int32_t, int64_t

// ----------------------------------------------------------------------
const char* pylith::meshio::PsetFileBinary::_HEADER = "pset unformatted";

// ----------------------------------------------------------------------
// Constructor with name of Pset file.
pylith::meshio::PsetFileBinary::PsetFileBinary(const char* filename,
                                               const bool flipEndian,
                                               const bool ioInt32,
                                               const bool isRecordHeader32Bit) :
    PsetFile(filename),
    _recordHeaderSize(isRecordHeader32Bit ? 4 : 8),
    _flipEndian(flipEndian),
    _ioInt32(ioInt32) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Default destructor
pylith::meshio::PsetFileBinary::~PsetFileBinary(void) { // destructor
} // destructor


// ----------------------------------------------------------------------
// Read binary Pset file.
void
pylith::meshio::PsetFileBinary::read(std::vector<Pset>* groups) { // read
    PYLITH_METHOD_BEGIN;

    assert(groups);

    pythia::journal::info_t info("meshio");

    std::ifstream fin(_filename.c_str(), std::ios::in);
    if (!(fin.is_open() && fin.good())) {
        std::ostringstream msg;
        msg
            << "Could not open binary Pset file '" << _filename
            << "' for reading.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Reading binary Pset file '" << _filename << "'." << pythia::journal::endl;

    _readHeader(fin);

    if (_ioInt32) {
        // Read number of psets
        int32_t numGroups = 0;
        fin.read((char*) &numGroups, sizeof(numGroups));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(numGroups));
        }
        assert(numGroups >= 0);
        groups->resize(numGroups);
        std::string extra = BinaryIO::readString(fin, 2*_recordHeaderSize);

        // Read groups
        info << pythia::journal::at(__HERE__)
             << "Reading " << numGroups << " point sets from file." << pythia::journal::endl;
        for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
            _readPset32(fin, &(*groups)[iGroup]);
        }
    } else {
        // Read number of psets
        int64_t numGroups = 0;
        fin.read((char*) &numGroups, sizeof(numGroups));
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(numGroups));
        }
        assert(numGroups >= 0);
        groups->resize(numGroups);
        std::string extra = BinaryIO::readString(fin, 2*_recordHeaderSize);

        // Read groups
        info << pythia::journal::at(__HERE__)
             << "Reading " << numGroups << " point sets from file." << pythia::journal::endl;
        for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
            _readPset64(fin, &(*groups)[iGroup]);
        }
    } // if/else

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write binary Pset file.
void
pylith::meshio::PsetFileBinary::write(const std::vector<Pset>& groups) { // write
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");

    std::ofstream fout(_filename.c_str(), std::ios::out);
    if (!(fout.is_open() && fout.good())) {
        std::ostringstream msg;
        msg
            << "Could not open binary Pset file '" << _filename
            << "' for writing.";
        throw std::runtime_error(msg.str());
    } // if

    info << pythia::journal::at(__HERE__)
         << "Writing binary Pset file '" << _filename << "'." << pythia::journal::endl;

    _writeHeader(fout);

    if (_ioInt32) {
        // Write number of groups
        int32_t numGroups = groups.size();
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(numGroups));
        }
        fout.write((char*) &numGroups, sizeof(numGroups));

        // Write groups
        info << pythia::journal::at(__HERE__)
             << "Writing " << numGroups << " point sets to file." << pythia::journal::endl;
        for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
            _writePset32(fout, groups[iGroup]);
        }
    } else {
        // Write number of groups
        int64_t numGroups = groups.size();
        if (_flipEndian) {
            BinaryIO::swapByteOrder((char*) &numGroups, 1, sizeof(numGroups));
        }
        fout.write((char*) &numGroups, sizeof(numGroups));

        // Write groups
        info << pythia::journal::at(__HERE__)
             << "Writing " << numGroups << " point sets to file." << pythia::journal::endl;
        for (int iGroup = 0; iGroup < numGroups; ++iGroup) {
            _writePset64(fout, groups[iGroup]);
        }
    } // if/else

    PYLITH_METHOD_END;
} // write


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_readHeader(std::ifstream& fin) { // _readHeader
    PYLITH_METHOD_BEGIN;

    std::string extra = BinaryIO::readString(fin, _recordHeaderSize);

    std::string header = BinaryIO::readString(fin, strlen(_HEADER));
    std::string headerE = _HEADER;
    headerE = headerE.substr(0, headerE.find_first_of(" "));
    if (headerE != header) {
        std::ostringstream msg;
        msg
            << "Header in binary Pset file '" << header
            << "' does not match anticipated header '" << headerE << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // _readHeader


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_writeHeader(std::ofstream& fout) { // _writeHeader
    PYLITH_METHOD_BEGIN;

    fout.write((char*) _HEADER, strlen(_HEADER));

    PYLITH_METHOD_END;
} // _writeHeader


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_readPset32(std::ifstream& fin,
                                            Pset* group) { // _readPset32
    PYLITH_METHOD_BEGIN;

    assert(group);

    pythia::journal::info_t info("meshio");

    group->name = BinaryIO::readString(fin, 32);

    int32_t id = 0;
    fin.read((char*) &id, sizeof(id));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));
    }

    int32_t size = 0;
    fin.read((char*) &size, sizeof(size));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &size, 1, sizeof(size));
    }
    assert(size >= 0);
    std::string extra = BinaryIO::readString(fin, 2*_recordHeaderSize);
    info << pythia::journal::at(__HERE__)
         << "Reading point set '" << group->name << "' with " << size
         << " points." << pythia::journal::endl;

    group->points.resize(size);
    fin.read((char*) &group->points[0], size*sizeof(int32_t));
    extra = BinaryIO::readString(fin, 2*_recordHeaderSize);
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &group->points[0], size, sizeof(int32_t));
    }

    group->points -= 1; // use zero base

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readPset32


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_writePset32(std::ofstream& fout,
                                             const Pset& group) { // _writePset32
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");
    const int size = group.points.size();
    info << pythia::journal::at(__HERE__)
         << "Writing point set '" << group.name << "' with " << size
         << " points." << pythia::journal::endl;

    fout.write((char*) group.name.c_str(), 32);

    int32_t id = group.id;
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));
    }
    fout.write((char*) &id, sizeof(id));

    int32_t sizeIO = size;
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &sizeIO, 1, sizeof(sizeIO));
    }
    fout.write((char*) &sizeIO, sizeof(sizeIO));

    std::valarray<int32_t> pointsIO(size);
    for (int i = 0; i < size; ++i) {
        pointsIO[i] = group.points[i];
    }
    pointsIO += 1; // switch from zero base to one base
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &pointsIO[0], size, sizeof(int32_t));
    }
    fout.write((char*) &pointsIO[0], size*sizeof(int32_t));

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _writePset32


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_readPset64(std::ifstream& fin,
                                            Pset* group) { // _readPset64
    PYLITH_METHOD_BEGIN;

    assert(group);

    pythia::journal::info_t info("meshio");

    group->name = BinaryIO::readString(fin, 32);

    int64_t id = 0;
    fin.read((char*) &id, sizeof(id));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));
    }

    int64_t size = 0;
    fin.read((char*) &size, sizeof(size));
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &size, 1, sizeof(size));
    }
    assert(size >= 0);
    std::string extra = BinaryIO::readString(fin, 2*_recordHeaderSize);
    info << pythia::journal::at(__HERE__)
         << "Reading point set '" << group->name << "' with " << size
         << " points." << pythia::journal::endl;

    group->points.resize(size);
    std::valarray<int64_t> pointsIO(size);
    fin.read((char*) &pointsIO[0], size*sizeof(int64_t));
    extra = BinaryIO::readString(fin, 2*_recordHeaderSize);
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &pointsIO[0], size, sizeof(int64_t));
    }

    for (int i = 0; i < size; ++i) {
        group->points[i] = pointsIO[i];
    }

    group->points -= 1; // use zero base

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _readPset64


// ----------------------------------------------------------------------
void
pylith::meshio::PsetFileBinary::_writePset64(std::ofstream& fout,
                                             const Pset& group) { // _writePset64
    PYLITH_METHOD_BEGIN;

    pythia::journal::info_t info("meshio");
    const int size = group.points.size();
    info << pythia::journal::at(__HERE__)
         << "Writing point set '" << group.name << "' with " << size
         << " points." << pythia::journal::endl;

    fout.write((char*) group.name.c_str(), 32);

    int64_t id = group.id;
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &id, 1, sizeof(id));
    }
    fout.write((char*) &id, sizeof(int));

    int64_t sizeIO = size;
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &sizeIO, 1, sizeof(sizeIO));
    }
    fout.write((char*) &sizeIO, sizeof(int));

    std::valarray<int64_t> pointsIO(size);
    pointsIO += 1; // switch from zero base to one base
    if (_flipEndian) {
        BinaryIO::swapByteOrder((char*) &pointsIO[0], size, sizeof(int64_t));
    }
    fout.write((char*) &pointsIO[0], size*sizeof(int64_t));

    info << pythia::journal::at(__HERE__)
         << "Done." << pythia::journal::endl;

    PYLITH_METHOD_END;
} // _writePset64


// End of file
