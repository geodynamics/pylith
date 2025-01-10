// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace meshio {
        class TestDataWriterVTK;
        class TestDataWriterVTK_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTK {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
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

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTK_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterVTK_Data(void);

    /// Destructor
    ~TestDataWriterVTK_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* timestepFilename; ///< Name of file with no data fields.
    const char* vertexFilename; ///< Name of file with vertex fields.
    const char* cellFilename; ///< Name of file with cell fields.

}; // TestDataWriterVTK_Data

// End of file
