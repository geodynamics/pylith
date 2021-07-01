// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestDataWriterVTK.hh
 *
 * @brief C++ TestDataWriterVTK object
 *
 * C++ unit testing for DataWriterVTK.
 */

#if !defined(pylith_meshio_testdatawritervtk_hh)
#define pylith_meshio_testdatawritervtk_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace meshio {
        class TestDataWriterVTK;

        class TestDataWriterVTK_Data;
    } // meshio
} // pylith

// =====================================================================================================================
class pylith::meshio::TestDataWriterVTK {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Check VTK file against archived file.
     *
     * @param filename Name of file to check.
     * @param t Time for file.
     * @param timeFormat Format of timestamp in filename.
     */
    static
    void checkFile(const char* filename,
                   const PylithScalar t,
                   const char* timeFormat);

}; // class TestDataWriterVTK

// =====================================================================================================================
class pylith::meshio::TestDataWriterVTK_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterVTK_Data(void);

    /// Destructor
    ~TestDataWriterVTK_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    const char* timestepFilename; ///< Name of file with no data fields.
    const char* vertexFilename; ///< Name of file with vertex fields.
    const char* cellFilename; ///< Name of file with cell fields.

}; // TestDataWriterVTK_Data

#endif // pylith_meshio_testdatawritervtk_hh

// End of file
