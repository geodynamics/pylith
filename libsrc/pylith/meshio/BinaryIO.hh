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

/**
 * @file libsrc/meshio/BinaryIO.hh
 *
 * @brief C++ object for general binary input/output operations.
 */

#if !defined(pylith_meshio_binaryio_hh)
#define pylith_meshio_binaryio_hh

#include "meshiofwd.hh" // forward declarations

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

#endif // pylith_meshio_binaryio_hh

// End of file
