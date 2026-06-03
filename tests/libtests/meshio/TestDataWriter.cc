// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestDataWriter.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::TestDataWriter_Data::TestDataWriter_Data(void) :
    meshFilename(nullptr),
    faultLabel(nullptr),
    faultId(-1),
    spaceDim(0),
    lengthScale(1.0),
    time(0.0),
    timeFormat(nullptr),
    fieldFn(nullptr) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::TestDataWriter_Data::~TestDataWriter_Data(void) {}


// End of file
