// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh" // forward declarations

#include <iosfwd>

class pylith::meshio::BinaryIO {
    // PUBLIC METHODS
    // //////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Read fixed length string from binary file.
     *
     * @param fin Input file stream
     * @param numChars Number of characters in string.
     */
    static
    std::string readString(std::ifstream& fin,
                           const int numChars);

    /** Change endian type by swapping byte order.
     *
     * @param vals Array of values
     * @param numVals Number of values
     * @param typesize Size of each value in bytes
     */
    static
    void swapByteOrder(char* vals,
                       const int numVals,
                       const int typesize);

}; // BinaryIO

// End of file
