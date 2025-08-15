// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestDataWriterVTK.hh" // Implementation of class methods

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <string.h> // USES strcmp()
#include <iostream> // USES std::cerr
#include <sstream> // USES std::ostringstream
#include <fstream> // USES std::ifstream

// ------------------------------------------------------------------------------------------------
// Check VTK file against archived file.
void
pylith::meshio::TestDataWriterVTK::checkFile(const char* filenameRoot,
                                             const PylithScalar t,
                                             const char* timeFormat) { // checkFile
    PYLITH_METHOD_BEGIN;

    const std::string& fileroot = filenameRoot;

    std::ostringstream buffer;
    const int indexExt = fileroot.find(".vtu");
    // Add time stamp to filename
    char sbuffer[256];
    snprintf(sbuffer, 256, timeFormat, t);
    std::string timestamp(sbuffer);
    const unsigned int pos = timestamp.find(".");
    if (pos != timestamp.length()) {
        timestamp.erase(pos, 1);
    }
    buffer << std::string(fileroot, 0, indexExt) << "_t" << timestamp << ".vtu";

    const std::string& filename = buffer.str();
    const std::string filenameE = "data/" + filename;

    std::ifstream fileInE(filenameE.c_str());
    if (!fileInE.is_open()) {
        std::cerr << "Could not open file '" << filenameE << "'." << std::endl;
    } // if
    assert(fileInE.is_open());

    std::ifstream fileIn(filename.c_str());
    if (!fileIn.is_open()) {
        std::cerr << "Could not open file '" << filename << "'." << std::endl;
    } // if
    assert(fileIn.is_open());

    const int maxLen = 256;
    char line[maxLen];
    char lineE[maxLen];
    const char* rawMarker = "AppendedData encoding=\"raw\"";
    int i = 1;
    while (!fileInE.eof()) {
        fileInE.getline(lineE, maxLen);
        fileIn.getline(line, maxLen);
        if (std::string(line).find(rawMarker) > 0) {
            break;
        } // if
        CHECK(std::string(lineE) == std::string(line));
        if (0 != strcmp(line, lineE)) {
            std::ostringstream msg;
            msg << "Mismatch in line " << i << " of file " << filename << ".\n"
                << "Expected: '" << lineE << "'\n"
                << "Actual: '" << line << "'";
            FAIL(msg.str());
        } // if
        ++i;
    } // while

    fileInE.close();
    fileIn.close();

    PYLITH_METHOD_END;
} // checkFile


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriterVTK_Data::TestDataWriterVTK_Data(void) :
    timestepFilename(NULL),
    vertexFilename(NULL),
    cellFilename(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriterVTK_Data::~TestDataWriterVTK_Data(void) {}


// End of file
