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

#include "BinaryIO.hh" // implementation of class methods

#include <fstream> // USES std::ifstream
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Read fixed length string from file.
std::string
pylith::meshio::BinaryIO::readString(std::ifstream& fin,
                                     const int numChars) { // readString
    std::string bstring = "";
    if (numChars > 0) {
        char* buffer = new char[numChars+1];
        buffer[0] = '\0';
        fin.read(buffer, sizeof(char)*numChars);
        buffer[numChars] = '\0';

        // get string from buffer
        std::string bufstring = buffer;
        delete[] buffer;buffer = 0;

        // remove whitespace
        const int iLast = bufstring.find_first_of(" ");
        bstring = bufstring.substr(0, iLast);
    } // if
    return std::string(bstring);
} // readString


// ---------------------------------------------------------------------------------------------------------------------
// Change endian type by swapping byte order.
void
pylith::meshio::BinaryIO::swapByteOrder(char* vals,
                                        const int numVals,
                                        const int typesize) { // swapByteOrder
    assert(vals);
    const int numSwaps = sizeof(typesize) / 2;
    for (int iVal = 0; iVal < numVals; ++iVal) {
        char* buf = (char*) (vals + iVal*typesize);
        for (int iSwap = 0, jSwap = typesize-1; iSwap < numSwaps; ++iSwap, --jSwap) {
            char tmp = buf[iSwap];
            buf[iSwap] = buf[jSwap];
            buf[jSwap] = tmp;
        } // for
    } // for
} // swapByteOrder


// End of file
